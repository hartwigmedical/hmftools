package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Integers.median;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.CN_SEGMENT_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

public class CopyNumberProfile
{
    private final PurityConfig mConfig;
    private final BufferedWriter mCnDataWriter;

    private final SampleData mSample;
    private final List<PurpleCopyNumber> mCopyNumbers;

    private final List<CopyNumberGcData> mCopyNumberGcRatios;

    public CopyNumberProfile(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mCnDataWriter = resultsWriter.getCnRatioWriter();
        mSample = sample;

        mCopyNumbers = Lists.newArrayList();

        mCopyNumberGcRatios = Lists.newArrayList();

        try
        {
            mCopyNumbers.addAll(PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(mConfig.getPurpleDir(sample.TumorId), mSample.TumorId)));
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
            return CnPurityResult.INVALID_RESULT;

        try
        {
            final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(mConfig.getCobaltDir(sampleId), sampleId);

            if(!Files.exists(Paths.get(cobaltFilename)))
            {
                CT_LOGGER.warn("sample({}) missing Cobalt ctDNA GC ratios file: {}", sampleId, cobaltFilename);
                return CnPurityResult.INVALID_RESULT;

            }

            Map<Chromosome,List<CobaltRatio>> cobaltRatios = CobaltRatioFile.readWithGender(cobaltFilename, null, true);

            buildCopyNumberGcRatios(cobaltRatios);

            if(mCnDataWriter != null)
            {
                mCopyNumberGcRatios.forEach(x -> writeCnSegmentData(mCnDataWriter, mConfig, mSample, sampleId, x));
            }

            double samplePloidy = purityContext.bestFit().ploidy();

            CnFitResult fitResult = CnPurityCalculator.calculatePurity(mCopyNumberGcRatios, samplePloidy);

            double fitPurityHigh = 0;
            double fitPurityLow = -1;
            // find a range by excluding each chromosome in turn
            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                List<CopyNumberGcData> excludedChrSegments = mCopyNumberGcRatios.stream()
                        .filter(x -> !RefGenomeFunctions.stripChrPrefix(x.Chromosome).equals(chromosome.toString()))
                        .collect(Collectors.toList());

                CnFitResult chrFitResult = CnPurityCalculator.calculatePurity(excludedChrSegments, samplePloidy);

                fitPurityLow = fitPurityLow < 0 ? chrFitResult.EstimatedPurity : min(fitPurityLow, chrFitResult.EstimatedPurity);
                fitPurityHigh = max(fitPurityHigh, chrFitResult.EstimatedPurity);
            }

            // if(!mPurityCalculator.valid())
            //    return INVALID_RESULT;

            CT_LOGGER.info(format("sample(%s) ploidy(%.4f) copy number segments(%d) estimated purity(%.6f)",
                    sampleId, samplePloidy, mCopyNumberGcRatios.size(), fitResult.EstimatedPurity));

            // calculate a median GC Ratio count and clonal percentage
            int totalRatios = 0;
            int clonalRatios = 0;

            List<Integer> segmentGcCounts = Lists.newArrayList();
            Map<Integer,Integer> cnLevelRatiosTotals = Maps.newHashMap();
            int maxCnLevelCount = 0;

            for(CopyNumberGcData cnSegment : mCopyNumberGcRatios)
            {
                int ratioCount = cnSegment.count();
                totalRatios += ratioCount;

                 if(cnSegment.IsValid)
                 {
                     clonalRatios += ratioCount;
                     segmentGcCounts.add(ratioCount);

                     Integer levelTotal = cnLevelRatiosTotals.get(cnSegment.CopyNumberLevel);
                     int newTotal = levelTotal != null ? levelTotal + ratioCount : ratioCount;
                     maxCnLevelCount = max(maxCnLevelCount, newTotal);
                     cnLevelRatiosTotals.put(cnSegment.CopyNumberLevel, newTotal);
                 }
            }

            double maxLevelPercent = clonalRatios > 0 ? maxCnLevelCount / (double)clonalRatios : 0;
            double anueploidyScore = 1 - maxLevelPercent;

            double clonalPercent = totalRatios > 0 ? clonalRatios / (double)totalRatios : 0;
            double medianGcRatioPerSegment = median(segmentGcCounts);

            return new CnPurityResult(
                    true, fitResult.FitCoefficient, fitResult.FitIntercept, fitResult.Residuals, fitResult.EstimatedPurity,
                    fitPurityLow, fitPurityHigh, anueploidyScore, clonalPercent,
                    mCopyNumberGcRatios.size(), mCopyNumberGcRatios.stream().mapToInt(x -> x.count()).sum(), medianGcRatioPerSegment);
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple and Cobalt copy-number data: {}", sampleId, e.toString());
            e.printStackTrace();
            return CnPurityResult.INVALID_RESULT;
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
        for(int i = 0; i <= PurityConstants.MAX_COPY_NUMBER; ++i)
        {
            if(abs(copyNumber - i) <= PurityConstants.CLONAL_COPY_NUMBER_MARGIN)
                return true;
        }

        return false;
    }

    private void plotCopyNumberGcRatioFit(final String sampleId)
    {
        plotCopyNumberGcRatioFit(mSample.PatientId, sampleId, mConfig);
    }

    public static BufferedWriter initialiseCnRatioWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(CN_SEGMENT_FILE_ID);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add("Chromosome").add("SegmentStart").add("SegmentEnd").add("CopyNumber");
            sj.add("GcRatioCount").add("GcRatioMedian").add("GcRatioMean");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise copy number segment file: {}", e.toString());
            return null;
        }
    }

    public static synchronized void writeCnSegmentData(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final CopyNumberGcData cnSegment)
    {
        if(writer == null)
            return;

        // writeGcRatioData(sampleId, cnSegment); // ratios are written before they are sorted for median calcs

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            addCommonFields(sj, config, sampleData, sampleId);

            sj.add(cnSegment.Chromosome);
            sj.add(String.valueOf(cnSegment.SegmentStart));
            sj.add(String.valueOf(cnSegment.SegmentEnd));
            sj.add(format("%.2f", cnSegment.CopyNumber));
            sj.add(String.valueOf(cnSegment.count()));
            sj.add(format("%.4f", cnSegment.median()));
            sj.add(format("%.4f", cnSegment.mean()));

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write copy number segment file: {}", e.toString());
        }
    }

    public static boolean plotCopyNumberGcRatioFit(final String patientId, final String sampleId, final PurityConfig config)
    {
        try
        {
            String summaryFile = config.formFilename(ResultsWriter.SUMMARY_FILE_ID);
            String cnSegmentsFile = config.formFilename(CN_SEGMENT_FILE_ID);

            if(!Files.exists(Paths.get(summaryFile)) || !Files.exists(Paths.get(cnSegmentsFile)))
            {
                CT_LOGGER.warn("plots missing required files: summary({}) segments({})", summaryFile, cnSegmentsFile);
                return false;
            }

            int runCode = RExecutor.executeFromClasspath(
                    "plots/CopyNumberGcRatioPlot.R", patientId, sampleId, summaryFile, cnSegmentsFile, config.PlotDir);

            return runCode == 0;
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to generate CN plot with R script: {}", e.toString());
            return false;
        }
    }
}
