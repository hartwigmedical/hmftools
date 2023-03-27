package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Integers.median;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.CnPurityResult.INVALID_RESULT;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.logging.log4j.Level;

public class CopyNumberProfile
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private PurityContext mPurityContext;
    private final List<PurpleCopyNumber> mCopyNumbers;

    private final List<CopyNumberGcData> mCopyNumberGcRatios;
    private final CnPurityCalculator mPurityCalculator;

    public CopyNumberProfile(final PurityConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;

        mCopyNumbers = Lists.newArrayList();
        mPurityContext = null;

        mCopyNumberGcRatios = Lists.newArrayList();
        mPurityCalculator = new CnPurityCalculator();

        try
        {
            mPurityContext = PurityContextFile.read(mConfig.PurpleDir, mConfig.TumorId);

            mCopyNumbers.addAll(PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDir, mConfig.TumorId)));
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple data: {}", mConfig.TumorId, e.toString());
        }
    }

    public List<CopyNumberGcData> copyNumberGcRatios() { return mCopyNumberGcRatios; }

    public CnPurityResult processSample(final String sampleId)
    {
        mCopyNumberGcRatios.clear();

        if(mPurityContext == null || mCopyNumbers.isEmpty())
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
                mCopyNumberGcRatios.forEach(x -> mResultsWriter.writeCnSegmentData(sampleId, x));

            double samplePloidy = mPurityContext.bestFit().ploidy();

            mPurityCalculator.calculatePurity(mCopyNumberGcRatios, samplePloidy);

            if(!mPurityCalculator.valid())
                return INVALID_RESULT;

            CT_LOGGER.info(format("sample(%s) ploidy(%.4f) copy number segments(%d) estimated purity(%.6f)",
                    sampleId, samplePloidy, mCopyNumberGcRatios.size(), mPurityCalculator.estimatedPurity()));

            double medianGcRatioPerSegment = calculateGcRatioCountMedian();

            return new CnPurityResult(
                    true, mPurityCalculator.fitCoefficient(), mPurityCalculator.fitIntercept(),
                    mPurityCalculator.residuals(), mPurityCalculator.estimatedPurity(),
                    mCopyNumberGcRatios.size(), mCopyNumberGcRatios.stream().mapToInt(x -> x.count()).sum(), medianGcRatioPerSegment);
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

            CopyNumberGcData cnSegment = new CopyNumberGcData(
                    copyNumber.chromosome(), copyNumber.start(), copyNumber.end(),
                    Doubles.round(copyNumber.averageTumorCopyNumber(), 2));

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

    private double calculateGcRatioCountMedian()
    {
        List<Integer> segmentGcCounts = Lists.newArrayList();
        mCopyNumberGcRatios.forEach(x -> segmentGcCounts.add(x.count()));
        return median(segmentGcCounts);
    }
}
