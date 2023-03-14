package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

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

public class CopyNumberProfile
{
    private final PurityConfig mConfig;
    private boolean mValid;

    public CopyNumberProfile(final PurityConfig config)
    {
        mConfig = config;
        mValid = true;
    }

    public boolean isValid() { return mValid; }

    public void loadSampleData(final String sampleId)
    {
        try
        {
            PurityContext purityContext = PurityContextFile.read(mConfig.PurpleDir, sampleId);

            List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDir, sampleId));

            final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(mConfig.CobaltDir, sampleId);

            Map<Chromosome,List<CobaltRatio>> cobaltRatios = CobaltRatioFile.readWithGender(cobaltFilename, null, true);

            // Multimap<Chromosome,GCProfile> gcProfile = GCProfileFactory.loadGCContent(mConfig.GcProfilePath);

            List<CopyNumberGcData> copyNumberGcRatios = buildCopyNumberGcRatios(cobaltRatios, copyNumbers);



        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple and Cobalt copy-number data: {}", sampleId, e.toString());
            e.printStackTrace();
            mValid = false;
        }
    }

    private List<CopyNumberGcData> buildCopyNumberGcRatios(
            final Map<Chromosome,List<CobaltRatio>> cobaltRatios, final List<PurpleCopyNumber> copyNumbers)
    {
        // expand the Purple copy numbers to segments to match GC profile
        List<CopyNumberGcData> copyNumberSegments = Lists.newArrayList();

        String currentChromosome = "";
        List<CobaltRatio> chrCobaltRatios = null;

        for(PurpleCopyNumber copyNumber : copyNumbers)
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
                    Doubles.round(copyNumber.averageTumorCopyNumber(), 1));

            segmentRatios.forEach(x -> cnSegment.addRatio(x.tumorGCRatio()));

            copyNumberSegments.add(cnSegment);

            CT_LOGGER.debug(format("segment(%s:%d - %d) copyNumber(%.4f) count(%d) mean(%.4f) median(%.4f)",
                    cnSegment.Chromosome, cnSegment.SegmentStart, cnSegment.SegmentEnd, cnSegment.CopyNumber,
                    cnSegment.count(), cnSegment.mean(), cnSegment.median()));
        }

        return copyNumberSegments;
    }
}
