package com.hartwig.hmftools.esvee.utils;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSvOrientation;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
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

import com.hartwig.hmftools.common.sv.gridss.GridssVcfTags;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
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
    private final Map<String, SAMRecord> mFilteredBamRecords = new HashMap<>();
    private final BufferedWriter mWriter;

    // For all gridss sam tags, see:
    // https://github.com/PapenfussLab/gridss/blob/master/src/main/java/au/edu/wehi/idsv/sam/SamTags.java
    private static final String SAM_TAG_READ_IDS = "ef";
    private static final String VCF_TAG_SV_TYPE = "EVENTTYPE";

    private static final String DELIMITER = TSV_DELIM;

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
        public final List<String> AssemblyIds;

        public Variant(final String chromosome, final int position, final byte orientation, String type, final List<String> assemblyIds)
        {
            Chromosome = chromosome;
            Position = position;
            Orientation = orientation;
            Type = type;
            AssemblyIds = assemblyIds;
        }

        public String toString()
        {
            return String.format("coord(%s:%s:%s) assemblyIds(%s)", Chromosome, Position, Orientation, AssemblyIds);
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

            Object assemblyIds = variantContext.getAttribute(GridssVcfTags.BEID);

            List<String> parsedAssemblyIds;
            if(assemblyIds == null)
                parsedAssemblyIds = new ArrayList<>();
            else if(assemblyIds instanceof String)
                parsedAssemblyIds = List.of((String) assemblyIds);
            else
                parsedAssemblyIds = (ArrayList<String>) assemblyIds;

            Variant variant = new Variant(
                    variantContext.getContig(),
                    variantContext.getStart(),
                    parseSvOrientation(variantContext),
                    (String) variantContext.getAttribute(VCF_TAG_SV_TYPE),
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

    private void filterBamAssemblies()
    {
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
                mFilteredBamRecords.put(bamAssemblyId, samRecord);
        }

        if(mFilteredBamRecords.size() == 0)
            SV_LOGGER.warn("BAM file({}) contained no assembly ids found in VCF file({})", mBamFile, mVcfFile);
        else
            SV_LOGGER.debug("Selected {} / {} assemblies with ids found in the VCF file", mFilteredBamRecords.size(), nTotalRecords);
    }

    private BufferedWriter initializeWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            String[] columnNames = { "VariantInfo", "VariantType", "ReadId", "AssemblyIds" };

            StringJoiner header = new StringJoiner(DELIMITER);
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

    private void writeVariantAssemblyReads(Variant variant)
    {
        try
        {
            String variantInfo = String.format("%s:%s:%s", variant.Chromosome, variant.Position, variant.Orientation);
            String assemblyIds = String.join(";", variant.AssemblyIds);

            // Need to dedup read ids as gridss sometimes reports a read id 2 (or more?) times per assembly
            Set<String> uniqueReadIds = new HashSet<>();

            for(String assemblyId : variant.AssemblyIds)
            {
                SAMRecord assembly = mFilteredBamRecords.get(assemblyId);

                String[] readIds = ((String) assembly.getAttribute(SAM_TAG_READ_IDS)).split(" ");
                uniqueReadIds.addAll(List.of(readIds));
            }

            for(String readId : uniqueReadIds)
            {
                String line = new StringJoiner(DELIMITER)
                        .add(variantInfo)
                        .add(variant.Type)
                        .add(readId)
                        .add(assemblyIds)
                        .toString();

                mWriter.write(line);
                mWriter.newLine();
            }
        }
        catch(Exception e)
        {
            SV_LOGGER.error("Failed to write assembly reads for variant({}): {}", variant.toString(), e.toString());
            System.exit(1);
        }
    }

    private void writeVariantAssemblyReads(int variantIndex)
    {
        writeVariantAssemblyReads(mVariants.get(variantIndex));
    }

    private void writeAssemblyReads()
    {
        SV_LOGGER.debug("Writing assembly reads to file: {}", mOutputFile);
        for(Variant variant : mVariants)
            writeVariantAssemblyReads(variant);
    }

    public void run()
    {
        try
        {
            SV_LOGGER.info("Extracting GRIDSS assembly reads from BAM file: {}", mBamFile);

            loadVcf();
            filterBamAssemblies();
            writeAssemblyReads();

            SV_LOGGER.info("Completed extracting assembly reads");
        }
        catch(Exception e)
        {
            SV_LOGGER.error("Failed to write assembly reads: {}", e.toString());
            System.exit(0);
        }
        finally
        {
            closeBufferedWriter(mWriter);
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
