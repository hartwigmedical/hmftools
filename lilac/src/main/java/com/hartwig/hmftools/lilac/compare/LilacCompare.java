package com.hartwig.hmftools.lilac.compare;

import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_ALIGN_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_INDELS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_FIT_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_HLA_Y;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_QC_STATUS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class LilacCompare
{
    private final String mSampleId;
    private final String mOutputDir;
    private final String mOutputId;
    private final String mOriginalDir;
    private final String mNewDir;

    private final double mPermittedDifference;
    private final double mPermittedDifferencePerc;

    private BufferedWriter mWriter;
    private int mMismatchCount;

    private static final String ORIGINAL_DIR = "original_dir";
    private static final String NEW_DIR = "new_dir";
    private static final String PERMITTED_DIFF = "diff_abs";
    private static final String PERMITTED_DIFF_PERC = "diff_perc";

    public LilacCompare(final ConfigBuilder configBuilder)
    {
        mOriginalDir = configBuilder.getValue(ORIGINAL_DIR);
        mNewDir = configBuilder.getValue(NEW_DIR);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);
        mSampleId = configBuilder.getValue(SAMPLE);

        mPermittedDifference = configBuilder.getDecimal(PERMITTED_DIFF);
        mPermittedDifferencePerc = configBuilder.getDecimal(PERMITTED_DIFF_PERC);

        mWriter = null;
        mMismatchCount = 0;
    }

    public void run()
    {
        LL_LOGGER.info("comparing Lilac output");
        LL_LOGGER.info("orig directory({})", mOriginalDir);
        LL_LOGGER.info("new directory({})", mNewDir);

        String outputFile = mOutputDir + mSampleId + LILAC_FILE_ID;

        if(mOutputId != null)
            outputFile += mOutputId + ".";

        outputFile += "compare.tsv";

        try
        {
            mWriter = createBufferedWriter(outputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("Field").add("OrigValue").add("NewValue");
            mWriter.write(sj.toString());
            mWriter.newLine();

            compareQcData();
            compareCandidateCoverage();
            compareFinalAlleles();

            mWriter.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write comparison data: {}", e.toString());
            System.exit(1);
        }

        LL_LOGGER.info("Lilac comparison complete, mismatches({})", mMismatchCount);
    }

    private void compareCandidateCoverage() throws IOException
    {
        // Score	HomozygousCount	CohortFrequency	RecoveryCount	WildcardCount	TotalCoverage
        // UniqueCoverage	SharedCoverage	WildCoverage	A1	A2	B1	B2	C1	C2
    }

    private void compareQcData() throws IOException
    {
        String origFile = LilacQcData.generateFilename(mOriginalDir, mSampleId);
        String newFile = LilacQcData.generateFilename(mNewDir, mSampleId);

        LilacQcData origQcData = LilacQcData.read(origFile);
        LilacQcData newQcData = LilacQcData.read(newFile);

        checkDifference(FLD_QC_STATUS, origQcData.status(), newQcData.status());
        checkDifference(FLD_TOTAL_FRAGS, origQcData.totalFragments(), newQcData.totalFragments());
        checkDifference(FLD_FIT_FRAGS, origQcData.fittedFragments(), newQcData.fittedFragments());
        checkDifference(FLD_DISC_ALIGN_FRAGS, origQcData.discardedAlignmentFragments(), newQcData.discardedAlignmentFragments());
        checkDifference(FLD_DISC_INDELS, origQcData.discardedIndels(), newQcData.discardedIndels());
        checkDifference(FLD_HLA_Y, origQcData.hlaYAllele(), newQcData.hlaYAllele());
    }

    private void compareFinalAlleles() throws IOException
    {
        // Allele	RefTotal	RefUnique	RefShared	RefWild	TumorTotal	TumorUnique	TumorShared	TumorWild	RnaTotal	RnaUnique	RnaShared	RnaWild
        // TumorCopyNumber	SomaticMissense	SomaticNonsenseOrFrameshiftSomaticSplice	SomaticSynonymous	SomaticInframeIndel
        String origFile = LilacAllele.generateFilename(mOriginalDir, mSampleId);
        String newFile = LilacAllele.generateFilename(mNewDir, mSampleId);

        List<LilacAllele> origAlleles = LilacAllele.read(origFile);
        List<LilacAllele> newAlleles = LilacAllele.read(newFile);

        if(origAlleles.size() != newAlleles.size())
        {
            writeMismatch("ALLELE_COUNTS", String.valueOf(origAlleles.size()), String.valueOf(newAlleles.size()));
            return;
        }

        List<LilacAllele> newAllelesToMatch = Lists.newArrayList(newAlleles);

        for(LilacAllele origAllele : origAlleles)
        {
            LilacAllele newAllele =
                    newAllelesToMatch.stream().filter(x -> x.allele().equals(origAllele.allele())).findFirst().orElse(null);

            if(newAllele != null)
            {
                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_REF_TOTAL), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_REF_UNIQUE), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_REF_SHARED), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_REF_WILD), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_TUMOR_TOTAL), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_SYNON), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_INDEL), origAllele.somaticMissense(), newAllele.somaticMissense());

                checkDifference(
                        format("%s_%s", origAllele.allele(), LilacAllele.FLD_NFS), origAllele.somaticMissense(), newAllele.somaticMissense());

                newAllelesToMatch.remove(newAllele);
            }
            else
            {
                writeMismatch("ALLELE", origAllele.allele(), "OTHER");
            }
        }

        for(LilacAllele newAllele : newAllelesToMatch)
        {
            writeMismatch("ALLELE", "OTHER", newAllele.allele());
        }
    }

    private void checkDifference(final String field, final String origValue, final String newValue) throws IOException
    {
        if(!origValue.equals(newValue))
            writeMismatch(field, origValue, newValue);
    }

    private void checkDifference(final String field, final int origValue, final int newValue) throws IOException
    {
        checkDifference(field, (double)origValue, (double)newValue);
    }

    private void checkDifference(final String field, final double origValue, final double newValue) throws IOException
    {
        if(origValue == newValue)
            return;

        double diff = abs(origValue - newValue);
        double diffPerc = diff / (double)max(origValue, newValue);

        if(diff > mPermittedDifference && diffPerc > mPermittedDifferencePerc)
        {
            writeMismatch(field, String.valueOf(origValue), String.valueOf(newValue));
        }
    }

    private void writeMismatch(final String field, final String origValue, final String newValue) throws IOException
    {
        ++mMismatchCount;
        mWriter.write(format("%s\t%s\t%s", field, origValue, newValue));
        mWriter.newLine();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(ORIGINAL_DIR, true, "Path to original Lilac output");
        configBuilder.addPath(NEW_DIR, true, "Path to new Lilac output");
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addDecimal(PERMITTED_DIFF, "Permitted numeric difference", 0.5);
        configBuilder.addDecimal(PERMITTED_DIFF_PERC, "Permitted numeric difference", 0.02);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        LilacCompare lilacCompare = new LilacCompare(configBuilder);
        lilacCompare.run();
    }
}
