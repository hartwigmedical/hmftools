package com.hartwig.hmftools.orange.algo.purple;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation.FULL_DEL;
import static com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation.FULL_GAIN;
import static com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation.PARTIAL_DEL;
import static com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation.PARTIAL_GAIN;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutableGermlineAmpDelFields;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

public final class GermlineGainDeletionFactory
{
    public static List<PurpleGainDeletion> createGermlineGainDeletions(
            final List<GermlineAmpDel> germlineAmpDels,
            final List<PurpleDriver> germlineDrivers,
            final List<GeneCopyNumber> geneCopyNumbers)
    {
        List<PurpleGainDeletion> gainDeletions = Lists.newArrayList();

        for(PurpleDriver driver : germlineDrivers)
        {
            if(driver.reportedStatus() != ReportedStatus.REPORTED)
                continue;

            if(driver.type() != PurpleDriverType.GERMLINE_DELETION && driver.type() != PurpleDriverType.GERMLINE_AMP)
                continue;

            List<GermlineAmpDel> matchedGermlineAmpDels = germlineAmpDels.stream()
                    .filter(x -> x.GeneName.equals(driver.gene())).collect(Collectors.toList());

            GeneCopyNumber geneCopyNumber = geneCopyNumbers.stream()
                    .filter(x -> x.GeneName.equals(driver.gene()) && driver.transcript().equals(x.TransName)).findFirst().orElse(null);

            if(matchedGermlineAmpDels.isEmpty() || geneCopyNumber == null)
                continue;

            gainDeletions.add(toGainDel(driver, matchedGermlineAmpDels, geneCopyNumber));
        }

        return gainDeletions;
    }

    private static PurpleGainDeletion toGainDel(
            final PurpleDriver driver, final List<GermlineAmpDel> germlineAmpDels, final GeneCopyNumber geneCopyNumber)
    {
        GermlineAmpDel firstGermlineAmpDel = germlineAmpDels.get(0);

        CopyNumberInterpretation interpretation;

        if(firstGermlineAmpDel.NormalStatus == AMPLIFICATION)
        {
            interpretation = firstGermlineAmpDel.IsPartial ? PARTIAL_GAIN : FULL_GAIN;
        }
        else
        {
            interpretation = firstGermlineAmpDel.IsPartial ? PARTIAL_DEL : FULL_DEL;
        }

        double minCopies = firstGermlineAmpDel.TumorCopyNumber;

        boolean coversTranscript = germlineAmpDels.stream().noneMatch(x -> x.IsPartial);

        double maxCopies = coversTranscript ?
                firstGermlineAmpDel.TumorCopyNumber : max(geneCopyNumber.maxCopyNumber(), firstGermlineAmpDel.TumorCopyNumber);

        double germlineMinCopies = firstGermlineAmpDel.GermlineCopyNumber;

        if(germlineAmpDels.size() > 1)
        {
            for(GermlineAmpDel germlineAmpDel : germlineAmpDels)
            {
                minCopies = min(germlineAmpDel.TumorCopyNumber, minCopies);
                maxCopies = max(germlineAmpDel.TumorCopyNumber, maxCopies);
                germlineMinCopies = min(germlineAmpDel.GermlineCopyNumber, germlineMinCopies);
            }
        }

        return ImmutablePurpleGainDeletion.builder()
                .driver(driver)
                .interpretation(interpretation)
                .chromosome(firstGermlineAmpDel.Chromosome)
                .chromosomeBand(firstGermlineAmpDel.ChromosomeBand)
                .germlineAmpDelFields(ImmutableGermlineAmpDelFields.builder()
                        .germlineStatus(PurpleConversion.convert(firstGermlineAmpDel.NormalStatus))
                        .somaticStatus(PurpleConversion.convert(firstGermlineAmpDel.TumorStatus))
                        .germlineMinCopyNumber(germlineMinCopies)
                        .build())
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .minMinorAlleleCopies(geneCopyNumber.MinMinorAlleleCopyNumber)
                .build();
    }
}
