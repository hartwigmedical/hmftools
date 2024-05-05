package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_FRAGMENTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_AA;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_FRAGS;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_CANDIDATE_NUC;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_READS;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_SOMATIC_VCF;
import static com.hartwig.hmftools.lilac.fragment.FragmentSource.REFERENCE;
import static com.hartwig.hmftools.lilac.fragment.FragmentSource.TUMOR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.coverage.HlaComplexFile;
import com.hartwig.hmftools.lilac.coverage.HlaYCoverage;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.FragmentSource;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
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
    private final BufferedWriter mFragmentWriter;
    private final BufferedWriter mReadWriter;
    private final Set<WriteType> mWriteTypes;
    private LilacVCF mVcfWriter;

    private static final String WRITE_TYPES = "write_types";
    private static final String WRITE_TYPES_ALL = "ALL";

    private enum WriteType
    {
        SUMMARY,
        FRAGMENTS,
        READS,
        REF_COUNTS;
    }

    public ResultsWriter(final LilacConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        mWriteTypes = Sets.newHashSet();

        String writeTypesStr = configBuilder.getValue(WRITE_TYPES);

        if(writeTypesStr.equals(WRITE_TYPES_ALL))
        {
            Arrays.stream(WriteType.values()).forEach(x -> mWriteTypes.add(x));
        }
        else
        {
            String[] types = writeTypesStr.split(",");
            Arrays.stream(types).forEach(x -> mWriteTypes.add(WriteType.valueOf(x)));
        }

        mFragmentWriter = mWriteTypes.contains(WriteType.FRAGMENTS) ? initialiseFragmentWriter() : null;
        mReadWriter = mWriteTypes.contains(WriteType.READS) ? initialiseReadWriter() : null;
        mVcfWriter = null;
    }

    public void close()
    {
        closeBufferedWriter(mFragmentWriter);
        closeBufferedWriter(mReadWriter);

        if(mVcfWriter != null)
            mVcfWriter.close();
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(WRITE_TYPES, false,
                "Write types: SUMMARY (default), FRAGMENTS, READS, REF_COUNTS or ALL", WriteType.SUMMARY.toString());
    }

    public void writeMainOutputs(
            final LilacQC summaryMetrics, final SolutionSummary solutionSummary, final List<ComplexCoverage> rankedComplexes)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        solutionSummary.write(LilacAllele.generateFilename(mConfig.OutputDir, mConfig.Sample));
        summaryMetrics.writefile(LilacQcData.generateFilename(mConfig.OutputDir, mConfig.Sample));

        HlaComplexFile.writeToFile(mConfig.formFileId(LILAC_FILE_CANDIDATE_COVERAGE), rankedComplexes);
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
            hlaYCoverage.writeAlleleCounts(mConfig.Sample);
        }
    }

    public void writeReferenceFragments(
            final List<ComplexCoverage> rankedComplexes, final List<Fragment> refNucleotideFrags, final List<FragmentAlleles> refFragAlleles)
    {
        if(!mWriteTypes.contains(WriteType.FRAGMENTS))
            return;

        HlaComplexFile.writeFragmentAssignment(mConfig.formFileId(LILAC_FILE_CANDIDATE_FRAGS), rankedComplexes, refFragAlleles);

        writeFragments(mConfig.tumorOnly() ? TUMOR : REFERENCE, refNucleotideFrags);
    }

    public void writeFailedSampleFileOutputs(final Map<String,int[]> geneBaseDepth, int medianBaseQuality)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        LL_LOGGER.info("writing failed-sample output to {}", mConfig.OutputDir);

        HaplotypeQC haplotypeQC = new HaplotypeQC(0, 0, 0, Lists.newArrayList());
        AminoAcidQC aminoAcidQC = new AminoAcidQC(0, 0);
        BamQC bamQC = new BamQC(0, 0, 0, geneBaseDepth);

        CoverageQC coverageQC = new CoverageQC(
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        LilacQC summaryMetrics = new LilacQC(
                0, "", 0, null,
                aminoAcidQC, bamQC, coverageQC, haplotypeQC, new SomaticVariantQC(0, 0));

        SolutionSummary solutionSummary = new SolutionSummary(
                null, null, null, null, null);

        solutionSummary.write(LilacAllele.generateFilename(mConfig.OutputDir, mConfig.Sample));
        summaryMetrics.writefile(LilacQcData.generateFilename(mConfig.OutputDir, mConfig.Sample));
    }

    private BufferedWriter initialiseFragmentWriter()
    {
        String filename = mConfig.formFileId(LILAC_FILE_FRAGMENTS);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);
            writer.write("Source\tReadId\tReadInfo\tGenes");
            writer.write("\tNucLociStart\tNucLociEnd\tAcidLociStart\tAcidLociEnd\tScope");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", filename, e.toString());
            return null;
        }
    }

    public void writeFragments(final FragmentSource source, final List<Fragment> fragments)
    {
        if(mFragmentWriter == null)
            return;

        try
        {
            for(Fragment fragment : fragments)
            {
                StringJoiner genesStr = new StringJoiner(ITEM_DELIM);
                fragment.getGenes().forEach(x -> genesStr.add(x));

                mFragmentWriter.write(String.format("%s\t%s\t%s\t%s",
                        source.toString(), fragment.id(), fragment.readInfo(), genesStr));

                mFragmentWriter.write(String.format("\t%d\t%d\t%d\t%d\t%s",
                        fragment.minNucleotideLocus(), fragment.maxNucleotideLocus(),
                        fragment.minAminoAcidLocus(), fragment.maxAminoAcidLocus(), fragment.scope()));

                mFragmentWriter.newLine();
            }
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write fragments: {}", e.toString());
        }
    }

    private BufferedWriter initialiseReadWriter()
    {
        String filename = mConfig.formFileId(LILAC_FILE_READS);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);
            writer.write("Source\tReadId\tReadInfo\tGenes");
            // writer.write("\tNucLociStart\tNucLociEnd\tAcidLociStart\tAcidLociEnd\tScope");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", filename, e.toString());
            return null;
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
