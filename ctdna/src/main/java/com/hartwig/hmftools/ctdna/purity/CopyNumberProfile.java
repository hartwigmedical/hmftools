package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Integers.median;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.CnPurityResult.INVALID_RESULT;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.CLONAL_COPY_NUMBER_MARGIN;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MAX_COPY_NUMBER;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.CN_SEGMENT_FILE_ID;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.SUMMARY_FILE_ID;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.logging.log4j.Level;

public class CopyNumberProfile
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final SampleData mSample;
    private final List<PurpleCopyNumber> mCopyNumbers;

    private final List<CopyNumberGcData> mCopyNumberGcRatios;
    private final CnPurityCalculator mPurityCalculator;

    public CopyNumberProfile(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;

        mCopyNumbers = Lists.newArrayList();

        mCopyNumberGcRatios = Lists.newArrayList();
        mPurityCalculator = new CnPurityCalculator();

        try
        {
            mCopyNumbers.addAll(PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDir, mSample.TumorId)));
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple data: {}", mSample.TumorId, e.toString());
        }
    }

    public CnPurityResult processSample(final String sampleId, final PurityContext purityContext)
    {
        mCopyNumberGcRatios.clear();

        if(purityContext == null || mCopyNumbers.isEmpty())
            return INVALID_RESULT;

        try
        {
            final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(mConfig.CobaltDir, sampleId);

            if(!Files.exists(Paths.get(cobaltFilename)))
            {
                CT_LOGGER.warn("sample({}) missing Cobalt ctDNA GC ratios file: {}", sampleId, cobaltFilename);
                return INVALID_RESULT;

            }

            Map<Chromosome,List<CobaltRatio>> cobaltRatios = CobaltRatioFile.readWithGender(cobaltFilename, null, true);

            buildCopyNumberGcRatios(cobaltRatios);

            if(mResultsWriter != null)
            {
                mCopyNumberGcRatios.forEach(x -> mResultsWriter.writeCnSegmentData(mSample.PatientId, sampleId, x));
            }

            double samplePloidy = purityContext.bestFit().ploidy();

            mPurityCalculator.calculatePurity(mCopyNumberGcRatios, samplePloidy);

            if(!mPurityCalculator.valid())
                return INVALID_RESULT;

            CT_LOGGER.info(format("sample(%s) ploidy(%.4f) copy number segments(%d) estimated purity(%.6f)",
                    sampleId, samplePloidy, mCopyNumberGcRatios.size(), mPurityCalculator.estimatedPurity()));

            // calculate a median GC Ratio count and clonal percentage
            int totalRatios = 0;
            int clonalRatios = 0;

            List<Integer> segmentGcCounts = Lists.newArrayList();
            Map<Integer,Integer> cnLevelRatioaTotals = Maps.newHashMap();
            int maxCnLevelCount = 0;

            for(CopyNumberGcData cnSegment : mCopyNumberGcRatios)
            {
                int ratioCount = cnSegment.count();
                totalRatios += ratioCount;

                 if(cnSegment.IsValid)
                 {
                     clonalRatios += ratioCount;
                     segmentGcCounts.add(ratioCount);

                     Integer levelTotal = cnLevelRatioaTotals.get(cnSegment.CopyNumberLevel);
                     int newTotal = levelTotal != null ? levelTotal + ratioCount : ratioCount;
                     maxCnLevelCount = max(maxCnLevelCount, newTotal);
                     cnLevelRatioaTotals.put(cnSegment.CopyNumberLevel, newTotal);
                 }
            }

            double maxLevelPercent = clonalRatios > 0 ? maxCnLevelCount / (double)clonalRatios : 0;
            double anueploidyScore = 1 - maxLevelPercent;

            double clonalPercent = totalRatios > 0 ? clonalRatios / (double)totalRatios : 0;
            double medianGcRatioPerSegment = median(segmentGcCounts);

            return new CnPurityResult(
                    true, mPurityCalculator.fitCoefficient(), mPurityCalculator.fitIntercept(),
                    mPurityCalculator.residuals(), mPurityCalculator.estimatedPurity(),
                    mCopyNumberGcRatios.size(), mCopyNumberGcRatios.stream().mapToInt(x -> x.count()).sum(), medianGcRatioPerSegment,
                    anueploidyScore, clonalPercent);
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple and Cobalt copy-number data: {}", sampleId, e.toString());
            e.printStackTrace();
            return INVALID_RESULT;
        }
    }

    private void buildCopyNumberGcRatios(final Map<Chromosome,List<CobaltRatio>> cobaltRatios)
    {
        // expand the Purple copy numbers to segments to match GC profile
        String currentChromosome = "";
        List<CobaltRatio> chrCobaltRatios = null;

        for(PurpleCopyNumber copyNumber : mCopyNumbers)
        {
            if(!currentChromosome.equals(copyNumber.chromosome()))
            {
                HumanChromosome chromosome = HumanChromosome.fromString(copyNumber.chromosome());

                if(chromosome.isAllosome())
                    continue;

                currentChromosome = copyNumber.chromosome();

                chrCobaltRatios = cobaltRatios.get(chromosome).stream().filter(x -> x.tumorGCRatio() >= 0).collect(Collectors.toList());
            }

            List<CobaltRatio> segmentRatios = chrCobaltRatios.stream()
                    .filter(x -> BaseRegion.positionWithin(x.position(), copyNumber.start(), copyNumber.end()))
                    .collect(Collectors.toList());

            if(segmentRatios.isEmpty())
                continue;

            boolean useCnSegment = useCopyNumberSegment(copyNumber.averageTumorCopyNumber());
            int cnLevel = (int)round(copyNumber.averageTumorCopyNumber());

            CopyNumberGcData cnSegment = new CopyNumberGcData(
                    copyNumber.chromosome(), copyNumber.start(), copyNumber.end(),
                    Doubles.round(copyNumber.averageTumorCopyNumber(), 2), useCnSegment, cnLevel);

            segmentRatios.forEach(x -> cnSegment.addRatio(new GcRatioData(x.position(), x.tumorGCRatio())));

            mCopyNumberGcRatios.add(cnSegment);

            if(CT_LOGGER.isTraceEnabled())
            {
                CT_LOGGER.trace(format("segment(%s:%d - %d) copyNumber(%.2f) count(%d) mean(%.4f) median(%.4f)",
                        cnSegment.Chromosome, cnSegment.SegmentStart, cnSegment.SegmentEnd, cnSegment.CopyNumber,
                        cnSegment.count(), cnSegment.mean(), cnSegment.median()));
            }
        }
    }

    private static boolean useCopyNumberSegment(double copyNumber)
    {
        for(int i = 0; i <= MAX_COPY_NUMBER; ++i)
        {
            if(abs(copyNumber - i) <= CLONAL_COPY_NUMBER_MARGIN)
                return true;
        }

        return false;
    }

    private void plotCopyNumberGcRatioFit(final String sampleId)
    {
        plotCopyNumberGcRatioFit(mSample.PatientId, sampleId, mConfig);
    }

    public static boolean plotCopyNumberGcRatioFit(final String patientId, final String sampleId, final PurityConfig config)
    {
        try
        {
            String summaryFile = config.formFilename(SUMMARY_FILE_ID);
            String cnSegmentsFile = config.formFilename(CN_SEGMENT_FILE_ID);

            int runCode = RExecutor.executeFromClasspath(
                    "plots/CopyNumberGcRatioPlot.R", patientId, sampleId, summaryFile, cnSegmentsFile, config.OutputDir);

            return runCode == 0;
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to generate CN plot with R script: {}", e.toString());
            return false;
        }
    }
}
