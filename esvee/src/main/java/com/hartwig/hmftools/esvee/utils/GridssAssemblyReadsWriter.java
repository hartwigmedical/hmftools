package com.hartwig.hmftools.esvee.utils;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.gridss.GridssVcfTags;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class GridssAssemblyReadsWriter
{
    public final String mBamFile;
    public final String mVcfFile;
    public final String mOutputFile;

    private final List<Variant> mVariants = new ArrayList<>();
    private final Map<String, SAMRecord> mAssemblies = new HashMap<>();
    private final BufferedWriter mWriter;
    int mReadsWritten = 0;

    // For all gridss sam tags, see:
    // https://github.com/PapenfussLab/gridss/blob/master/src/main/java/au/edu/wehi/idsv/sam/SamTags.java
    private static final String SAM_TAG_READ_IDS = "ef";

    private static final String VCF_TAG_SV_TYPE = "EVENTTYPE";
    private static final String VCF_TAG_SV_ID = "EVENT";
    private static final String VCF_TAG_CIPOS = "CIPOS";

    private static final String BAM_FILE = "bam_file";
    private static final String VCF_FILE = "vcf_file";
    private static final String OUTPUT_TSV = "output_tsv";

    public GridssAssemblyReadsWriter(String bamFile, String vcfFile, String outputFile)
    {
        mBamFile = bamFile;
        mVcfFile = vcfFile;
        mOutputFile = outputFile;
        mWriter = initializeWriter();
    }

    private static class Variant
    {
        public final String Chromosome;
        public final int Position;
        public final byte Orientation;
        public final String Type;
        public final int CiLower;
        public final int CiUpper;
        public final String SvId;
        public final String BreakendId;
        public final List<String> AssemblyIds;

        public String IntersectingAssemblyId;

        public Variant(final String chromosome, final int position, final byte orientation, String type,
                final int ciLower, final int ciUpper,
                String svId, String breakendId,
                final List<String> assemblyIds)
        {
            Chromosome = chromosome;
            Position = position;
            Orientation = orientation;
            Type = type;
            CiLower = ciLower;
            CiUpper = ciUpper;
            SvId = svId;
            BreakendId = breakendId;
            AssemblyIds = assemblyIds;

            IntersectingAssemblyId = null;
        }

        public String toString()
        {
            return String.format("%s:%s:%s", Chromosome, Position, Orientation);
        }
    }

    private void loadVcf()
    {
        VcfFileReader reader = new VcfFileReader(mVcfFile);

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        int nTotalVariants = 0;
        int nSelectedVariants = 0;
        for(VariantContext variantContext : reader.iterator())
        {
            nTotalVariants++;
            if(!filter.test(variantContext))
                continue;

            //
            List<String> cipos = (List<String>) variantContext.getAttribute(VCF_TAG_CIPOS);
            int ciLower = 0;
            int ciUpper = 0;
            if(cipos != null)
            {
                ciLower = Integer.parseInt( cipos.get(0) );
                ciUpper = Integer.parseInt( cipos.get(1) );
            }

            //
            Object assemblyIds = variantContext.getAttribute(GridssVcfTags.BEID);
            List<String> parsedAssemblyIds;
            if(assemblyIds == null)
                parsedAssemblyIds = new ArrayList<>();
            else if(assemblyIds instanceof String)
                parsedAssemblyIds = List.of((String) assemblyIds);
            else
                parsedAssemblyIds = (List<String>) assemblyIds;

            byte orientation = StructuralVariantFactory.parseSvOrientation(variantContext);
            if(orientation == 0)
                orientation = StructuralVariantFactory.parseSingleOrientation(variantContext);

            Variant variant = new Variant(
                    variantContext.getContig(),
                    variantContext.getStart(),
                    orientation,
                    (String) variantContext.getAttribute(VCF_TAG_SV_TYPE),
                    ciLower,
                    ciUpper,
                    (String) variantContext.getAttribute(VCF_TAG_SV_ID),
                    variantContext.getID(),
                    parsedAssemblyIds
            );

            mVariants.add(variant);
            nSelectedVariants++;
        }

        if(mVariants.size() == 0)
            SV_LOGGER.warn("No variants found in VCF file: {}", mVcfFile);
        else
            SV_LOGGER.debug("Loaded {} / {} variants from VCF file: {}", nSelectedVariants, nTotalVariants, mVcfFile);
    }

    private void loadRelevantAssemblies()
    {
        // Assemblies (i.e. BAM records) cannot be retrieved by id which are contained in the QNAME (first) field.
        // Instead, we need to iterate over the BAM records and store those with assembly IDs that are found in the VCF.

        if(mVariants.size() == 0)
            return;

        HashSet<String> vcfAssemblyIds = new HashSet<>();
        for(Variant variant : mVariants)
        {
            if(variant.AssemblyIds.size() == 0)
                continue;

            vcfAssemblyIds.addAll(variant.AssemblyIds);
        }

        if(vcfAssemblyIds.size() == 0)
        {
            SV_LOGGER.warn("No assembly ids found in vcfFile: {}", mVcfFile);
            return;
        }

        SamReader samReader = SamReaderFactory
                .makeDefault()
                .validationStringency(ValidationStringency.STRICT)
                .open(new File(mBamFile));

        int nTotalRecords = 0;
        for(SAMRecord samRecord : samReader)
        {
            nTotalRecords++;

            String bamAssemblyId = samRecord.getReadName();
            if(vcfAssemblyIds.contains(bamAssemblyId))
                mAssemblies.put(bamAssemblyId, samRecord);
        }

        if(mAssemblies.size() == 0)
            SV_LOGGER.warn("BAM file({}) contained no assembly ids found in VCF file({})", mBamFile, mVcfFile);
        else
            SV_LOGGER.debug("Selected {} / {} assemblies with ids found in the VCF file", mAssemblies.size(), nTotalRecords);
    }

    private void findIntersectingAssemblies()
    {
        for(Variant variant : mVariants){

            String intersectingAssemblyId = null;

            for(String assemblyId : variant.AssemblyIds)
            {
                SAMRecord assembly = mAssemblies.get(assemblyId);

                boolean assemblyIntersectsBreakend = assembly.getReferenceName().equals(variant.Chromosome) &&
                        assembly.getAlignmentStart() <= variant.Position &&
                        assembly.getAlignmentEnd() >= variant.Position;

                if(assemblyIntersectsBreakend)
                {
                    intersectingAssemblyId = assemblyId;
                    break;
                }
            }

            variant.IntersectingAssemblyId = intersectingAssemblyId;
        }
    }

    private BufferedWriter initializeWriter()
    {
        try
        {
            BufferedWriter writer = FileWriterUtils.createBufferedWriter(mOutputFile, false);

            String[] columnNames = {
                    "Chromosome", "Position", "Orientation", "VariantType",
                    "CiLower", "CiUpper",
                    "SvId", "BreakendId",
                    "ReadId", "AssemblyId", "IsIntersectingAssembly"
            };

            StringJoiner header = new StringJoiner(TSV_DELIM);
            for(String columnName : columnNames)
                header.add(columnName);

            writer.write(header.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("Failed to initialise assembly reads writer: {}", e.toString());
            return null;
        }
    }

    private void writeVariantAssemblyReads(Variant variant, String assemblyId, boolean isIntersectingAssembly)
    {
        try
        {
            SAMRecord assembly = mAssemblies.get(assemblyId);

            // Need to dedup read ids as gridss sometimes reports a read id 2 (or more?) times per assembly
            String[] readIds = String.valueOf(assembly.getAttribute(SAM_TAG_READ_IDS)).split(" ");
            Set<String> uniqueReadIds = new HashSet<>(List.of(readIds));

            for(String readId : uniqueReadIds)
            {
                String line = String.join(
                        TSV_DELIM,
                        variant.Chromosome,
                        String.valueOf(variant.Position),
                        String.valueOf(variant.Orientation),
                        variant.Type,
                        String.valueOf(variant.CiLower),
                        String.valueOf(variant.CiUpper),
                        variant.SvId,
                        variant.BreakendId,
                        readId,
                        assembly.getReadName(),
                        String.valueOf(isIntersectingAssembly)
                );

                mWriter.write(line);
                mWriter.newLine();

                mReadsWritten++;
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("Failed to write assembly reads for variant({})", variant);
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void writeAssemblyReads()
    {
        SV_LOGGER.debug("Writing assembly reads to file: {}", mOutputFile);

        int variantsIntersectingAssembly = 0;

        for(Variant variant : mVariants)
        {
            if(variant.IntersectingAssemblyId != null)
            {
                writeVariantAssemblyReads(variant, variant.IntersectingAssemblyId, true);
                variantsIntersectingAssembly++;
            }
            else
            {
                SV_LOGGER.warn("variant({}) type({}) had no intersecting assemblies({}). Writing reads for all assemblies", variant, variant.Type, String.join(", ", variant.AssemblyIds));
                for(String assemblyId : variant.AssemblyIds)
                    writeVariantAssemblyReads(variant, assemblyId, false);
            }
        }

        SV_LOGGER.info("Found {} / {} variants intersecting one or more assemblies", variantsIntersectingAssembly, mVariants.size());
        SV_LOGGER.info("Wrote {} reads to file: {}", mReadsWritten, mOutputFile);
    }

    public void run()
    {
        try
        {
            SV_LOGGER.info("Extracting GRIDSS assembly reads from BAM file: {}", mBamFile);

            loadVcf();
            loadRelevantAssemblies();
            findIntersectingAssemblies();
            writeAssemblyReads();

            SV_LOGGER.info("Completed extracting assembly reads");
        }
        catch(Exception e)
        {
            SV_LOGGER.error("Failed to write assembly reads:");
            e.printStackTrace();
            System.exit(1);
        }
        finally
        {
            FileWriterUtils.closeBufferedWriter(mWriter);
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(VCF_FILE, true, "Input GRIDSS or PURPLE VCF file");
        configBuilder.addPath(BAM_FILE, true, "Input GRIDSS assemblies BAM file");
        configBuilder.addConfigItem(OUTPUT_TSV, true, "Output TSV file");
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GridssAssemblyReadsWriter writer = new GridssAssemblyReadsWriter(
                configBuilder.getValue(BAM_FILE),
                configBuilder.getValue(VCF_FILE),
                configBuilder.getValue(OUTPUT_TSV)
        );

        writer.run();
    }
}
