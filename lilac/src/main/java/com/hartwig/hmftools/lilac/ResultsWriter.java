package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.LilacConstants.CURRENT_GENES;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_AA;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_FRAGS;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_NUC;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_FRAGMENTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_SOMATIC_VCF;
import static com.hartwig.hmftools.lilac.fragment.FragmentSource.REFERENCE;
import static com.hartwig.hmftools.lilac.fragment.FragmentSource.TUMOR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;
import com.hartwig.hmftools.common.utils.file.FileLock;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.coverage.HlaComplexFile;
import com.hartwig.hmftools.lilac.coverage.HlaYCoverage;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.FragmentSource;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.qc.AminoAcidQC;
import com.hartwig.hmftools.lilac.qc.BamQC;
import com.hartwig.hmftools.lilac.qc.CoverageQC;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.qc.LilacQC;
import com.hartwig.hmftools.lilac.qc.SolutionSummary;
import com.hartwig.hmftools.lilac.qc.SomaticVariantQC;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.variant.LilacVCF;

import htsjdk.variant.variantcontext.VariantContext;

public class ResultsWriter
{
    private final LilacConfig mConfig;
    private final FileLock mFragmentWriter;
    private final Set<WriteType> mWriteTypes;
    private LilacVCF mVcfWriter;

    public static final String WRITE_TYPES = "write_types";
    private static final String WRITE_TYPES_ALL = "ALL";

    private enum WriteType
    {
        SUMMARY,
        FRAGMENTS,
        REF_COUNTS
    }

    public ResultsWriter(final LilacConfig config, final String writeTypesStr)
    {
        mConfig = config;
        mWriteTypes = Sets.newHashSet();
        if(writeTypesStr.equals(WRITE_TYPES_ALL))
        {
            mWriteTypes.addAll(Arrays.asList(WriteType.values()));
        }
        else
        {
            String[] types = writeTypesStr.split(",");
            Arrays.stream(types).forEach(x -> mWriteTypes.add(WriteType.valueOf(x)));
        }

        mFragmentWriter = mWriteTypes.contains(WriteType.FRAGMENTS) ? initialiseFragmentWriter() : null;
        mVcfWriter = null;
    }

    public void close()
    {
        try
        {
            if(mFragmentWriter != null)
                mFragmentWriter.close();

            if(mVcfWriter != null)
                mVcfWriter.close();
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(WRITE_TYPES, false,
                "Write types: SUMMARY (default), FRAGMENTS, READS, REF_COUNTS or ALL", WriteType.SUMMARY.toString());
    }

    public void writeMainOutputs(
            final LilacQC summaryMetrics, final SolutionSummary solutionSummary, final Iterable<ComplexCoverage> rankedComplexes)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        solutionSummary.write(LilacAllele.generateFilename(mConfig.OutputDir, mConfig.Sample));
        summaryMetrics.writefile(LilacQcData.generateFilename(mConfig.OutputDir, mConfig.Sample));

        HlaComplexFile.writeToFile(mConfig.formFileId(LILAC_FILE_CANDIDATE_COVERAGE, false), CURRENT_GENES, rankedComplexes);
    }

    public void writeDetailedOutputs(
            final SequenceCount refAminoAcidCounts, final SequenceCount refNucleotideCounts,
            final AminoAcidFragmentPipeline aminoAcidPipeline, final HlaYCoverage hlaYCoverage)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        if(mWriteTypes.contains(WriteType.REF_COUNTS))
        {
            refAminoAcidCounts.writeVertically(mConfig.formFileId(LILAC_FILE_CANDIDATE_AA));
            refNucleotideCounts.writeVertically(mConfig.formFileId(LILAC_FILE_CANDIDATE_NUC));
            aminoAcidPipeline.writeCounts(mConfig);
            if(hlaYCoverage != null)
                hlaYCoverage.writeAlleleCounts(mConfig.Sample);
        }
    }

    public void writeReferenceFragments(
            final Iterable<ComplexCoverage> rankedComplexes, final Iterable<Fragment> refNucleotideFrags,
            final List<FragmentAlleles> refFragAlleles)
    {
        if(!mWriteTypes.contains(WriteType.FRAGMENTS))
            return;

        String filename = mConfig.formFileId(LILAC_FILE_CANDIDATE_FRAGS, false);
        try(FileLock file = FileLock.create(new File(filename)))
        {
            List<String> existingLines = Lists.newArrayList();
            BufferedReader reader = file.getBufferedReader();
            reader.readLine();
            String line = reader.readLine();
            while(line != null)
            {
                String geneStr = line.split("\\*")[0];
                HlaGene gene = HlaGene.fromString(geneStr);
                if(gene != null && !CURRENT_GENES.contains(gene))
                    existingLines.add(line);

                line = reader.readLine();
            }

            file.clear();
            BufferedWriter writer = file.getBufferedWriter();
            HlaComplexFile.writeFragmentAssignment(writer, rankedComplexes, refFragAlleles);
            for(String existingLine : existingLines)
            {
                writer.write(existingLine);
                writer.newLine();
            }

            writer.flush();
        }
        catch(Exception e)
        {
            LL_LOGGER.error("failed to update {}: {}", filename, e.toString());
        }

        writeFragments(mConfig.tumorOnly() ? TUMOR : REFERENCE, refNucleotideFrags);
    }

    public void writeFailedSampleFileOutputs(final Map<HlaGene, int[]> geneBaseDepth)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        LL_LOGGER.info("writing failed-sample output to {}", mConfig.OutputDir);

        HaplotypeQC haplotypeQC = new HaplotypeQC(0, 0, 0, Lists.newArrayList());
        AminoAcidQC aminoAcidQC = new AminoAcidQC(0, 0);
        BamQC bamQC = new BamQC(0, 0, 0, geneBaseDepth);

        Map<HlaGene, Integer> countsByGene = Maps.newHashMap();
        for(HlaGene gene : geneBaseDepth.keySet())
            countsByGene.put(gene, 0);

        CoverageQC coverageQC = new CoverageQC(countsByGene, 0, 0, 0, 0, 0, 0, 0);

        LilacQC summaryMetrics = new LilacQC(
                CURRENT_GENES, 0, "", null, null,
                aminoAcidQC, bamQC, coverageQC, haplotypeQC, new SomaticVariantQC(0, 0));

        SolutionSummary solutionSummary = new SolutionSummary(
                CURRENT_GENES, null, null, null, null, null);

        solutionSummary.write(LilacAllele.generateFilename(mConfig.OutputDir, mConfig.Sample));
        summaryMetrics.writefile(LilacQcData.generateFilename(mConfig.OutputDir, mConfig.Sample));
    }

    private FileLock initialiseFragmentWriter()
    {
        String filename = mConfig.formFileId(LILAC_FILE_FRAGMENTS, false);
        FileLock fileLock = FileLock.create(new File(filename));
        if(fileLock == null)
        {
            LL_LOGGER.error("failed to update {}", filename);
            return null;
        }

        List<String> existingRecords = existingFragments(fileLock);
        try
        {
            fileLock.clear();
            BufferedWriter writer = fileLock.getBufferedWriter();
            writer.write("Source\tReadId\tReadInfo\tGenes");
            writer.write("\tNucLociStart\tNucLociEnd\tAcidLociStart\tAcidLociEnd\tScope");
            writer.newLine();
            for(String existingRecord : existingRecords)
            {
                writer.write(existingRecord);
                writer.newLine();
            }

            writer.flush();
            return fileLock;
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to update {}: {}", filename, e.toString());
            return null;
        }
    }

    private static List<String> existingFragments(final FileLock file)
    {
        List<String> existingFragmentRecords = Lists.newArrayList();
        try
        {
            BufferedReader reader = file.getBufferedReader();
            reader.readLine();
            String line = reader.readLine();
            while(line != null)
            {
                String[] fields = line.split("\t");
                if(fields.length < 4)
                {
                    LL_LOGGER.error("existing data in lilac.{} file is malformed", LILAC_FILE_FRAGMENTS);
                    return Collections.emptyList();
                }

                String[] geneStrs = fields[3].split(ITEM_DELIM);
                HlaGene firstGene = HlaGene.fromString(geneStrs[0]);
                if(firstGene != null && !CURRENT_GENES.contains(firstGene))
                    existingFragmentRecords.add(line);

                line = reader.readLine();
            }

            return existingFragmentRecords;
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read existing data in lilac.{} file: {}", LILAC_FILE_FRAGMENTS, e.toString());
            return Collections.emptyList();
        }
    }

    public void writeFragments(final FragmentSource source, final Iterable<Fragment> fragments)
    {
        if(mFragmentWriter == null)
            return;

        try
        {
            BufferedWriter writer = mFragmentWriter.getBufferedWriter();
            for(Fragment fragment : fragments)
            {
                StringJoiner genesStr = new StringJoiner(ITEM_DELIM);
                fragment.genes().forEach(x -> genesStr.add(x.toString()));

                writer.write(String.format("%s\t%s\t%s\t%s",
                        source.toString(), fragment.id(), fragment.readInfo(), genesStr));

                writer.write(String.format("\t%d\t%d\t%d\t%d\t%s",
                        fragment.minNucleotideLocus(), fragment.maxNucleotideLocus(),
                        fragment.minAminoAcidLocus(), fragment.maxAminoAcidLocus(), fragment.scope()));

                writer.newLine();
            }

            writer.flush();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write fragments: {}", e.toString());
        }
    }

    public void writeVariant(final VariantContext context, final List<HlaAllele> alleles)
    {
        if(mConfig.OutputDir.isEmpty() || context != null)
            return;

        if(mVcfWriter == null)
        {
            final VersionInfo version = fromAppName(APP_NAME);

            mVcfWriter = new LilacVCF(
                    mConfig.formFileId(LILAC_FILE_SOMATIC_VCF),
                    mConfig.SomaticVariantsFile).writeHeader(version.version());
        }

        mVcfWriter.writeVariant(context, alleles);
    }
}
