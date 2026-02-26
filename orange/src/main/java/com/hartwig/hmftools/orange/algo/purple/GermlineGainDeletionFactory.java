package com.hartwig.hmftools.orange.algo.purple;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.ImmutableGermlineAmpDelFields;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.Nullable;

public final class GermlineGainDeletionFactory
{
    public static List<PurpleGainDeletion> createGermlineGainDeletions(
            final List<GermlineAmpDel> germlineAmpDels, final List<PurpleDriver> germlineDrivers,
            final List<GeneCopyNumber> geneCopyNumbers, @Nullable final IsofoxData isofoxData)
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

            GeneExpression geneExpression = isofoxData != null ? isofoxData.geneExpressions().stream()
                    .filter(x -> x.geneName().equals(driver.gene())).findFirst().orElse(null) : null;

            gainDeletions.add(toGainDel(driver, matchedGermlineAmpDels, geneCopyNumber, geneExpression));
        }

        return gainDeletions;
    }

    private static PurpleGainDeletion toGainDel(
            final PurpleDriver driver, final List<GermlineAmpDel> germlineAmpDels, final GeneCopyNumber geneCopyNumber,
            @Nullable final GeneExpression geneExpression)
    {
        GermlineAmpDel firstGermlineAmpDel = germlineAmpDels.get(0);

        String exonRange = firstGermlineAmpDel.IsPartial ?
                format("%d_%d", firstGermlineAmpDel.ExonStart, firstGermlineAmpDel.ExonEnd) : "FULL";

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

        Double tpm = null;
        Double tpmPercentile = null;
        Double tpmFoldChange = null;

        if(geneExpression != null)
        {
            tpm = geneExpression.tpm();
            tpmPercentile = geneExpression.percentileCohort();
            tpmFoldChange = geneExpression.medianTpmCohort() > 0 ? tpm / geneExpression.medianTpmCohort() : 0;
        }

        return ImmutablePurpleGainDeletion.builder()
                .driver(driver)
                .chromosome(firstGermlineAmpDel.Chromosome)
                .chromosomeBand(firstGermlineAmpDel.ChromosomeBand)
                .germlineAmpDelFields(ImmutableGermlineAmpDelFields.builder()
                        .germlineStatus(PurpleConversion.convert(firstGermlineAmpDel.NormalStatus))
                        .somaticStatus(PurpleConversion.convert(firstGermlineAmpDel.TumorStatus))
                        .germlineMinCopyNumber(germlineMinCopies)
                        .build())
                .minCopyNumber(minCopies)
                .maxCopyNumber(maxCopies)
                .relativeCopyNumber(geneCopyNumber.RelativeMinCopyNumber)
                .exonRange(exonRange)
                .minMinorAlleleCopies(geneCopyNumber.MinMinorAlleleCopyNumber)
                .tpm(tpm)
                .tpmPercentile(tpmPercentile)
                .tpmFoldChange(tpmFoldChange)
                .build();
    }
}
