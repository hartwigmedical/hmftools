package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutableGermlineAmpDelFields;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

public class GermlineGainDeletionFactory
{
    private final EnsemblDataCache mEnsemblDataCache;

    public GermlineGainDeletionFactory(final EnsemblDataCache ensemblDataCache)
    {
        mEnsemblDataCache = ensemblDataCache;
    }

    public List<PurpleGainDeletion> createGermlineGainDeletions(
            final List<GermlineAmpDel> germlineAmpDels,
            final List<PurpleDriver> germlineDrivers,
            final List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        // get all the drivers
        List<PurpleDriver> reportableGermlineAmpDelDrivers = germlineDrivers.stream()
                .filter(o -> o.reportedStatus() == ReportedStatus.REPORTED)
                .filter(o -> o.type() == PurpleDriverType.GERMLINE_DELETION || o.type() == PurpleDriverType.GERMLINE_AMP)
                .toList();

        return reportableGermlineAmpDelDrivers.stream()
                .map(driver -> {
                    List<GermlineAmpDel> ampDels =
                            germlineAmpDels.stream().filter(o -> o.GeneName.equals(driver.gene())).toList();
                    GeneCopyNumber somaticGeneCopyNumber = GermlineGainDeletionUtil.findGeneCopyNumberForGeneTranscript(
                            driver.gene(), driver.transcript(), allSomaticGeneCopyNumbers);
                    return toGainDel(driver, ampDels, somaticGeneCopyNumber);
                })
                .toList();
    }

    private PurpleGainDeletion toGainDel(
            final PurpleDriver driver, final List<GermlineAmpDel> deletionsForGene, final GeneCopyNumber somaticGeneCopyNumber)
    {
        TranscriptData canonicalTranscript = GermlineGainDeletionUtil.findTranscript(driver.gene(), driver.transcript(), mEnsemblDataCache);
        CopyNumberInterpretation interpretation = GermlineGainDeletionUtil.deletionsCoverTranscript(deletionsForGene, canonicalTranscript)
                ? CopyNumberInterpretation.FULL_DEL
                : CopyNumberInterpretation.PARTIAL_DEL;
        double minCopies = GermlineGainDeletionUtil.getSomaticMinCopyNumber(deletionsForGene);
        double maxCopies = GermlineGainDeletionUtil.getSomaticMaxCopyNumber(deletionsForGene, somaticGeneCopyNumber, canonicalTranscript);
        double germlineMinCopies = GermlineGainDeletionUtil.getGermlineMinCopyNumber(deletionsForGene);

        GermlineStatus germlineStatus = GermlineGainDeletionUtil.getGermlineGainDelStatus(deletionsForGene);
        GermlineStatus somaticStatus = GermlineGainDeletionUtil.getSomaticGainDelStatus(deletionsForGene);

        String chromosome = GermlineGainDeletionUtil.getChromosome(deletionsForGene);
        String chromosomeBand = GermlineGainDeletionUtil.getChromosomeBand(deletionsForGene);

        return ImmutablePurpleGainDeletion.builder()
                .driver(driver)
                .interpretation(interpretation)
                .chromosome(chromosome)
                .chromosomeBand(chromosomeBand)
                .germlineAmpDelFields(ImmutableGermlineAmpDelFields.builder()
                        .germlineStatus(PurpleConversion.convert(germlineStatus))
                        .somaticStatus(PurpleConversion.convert(somaticStatus))
                        .germlineMinCopyNumber(germlineMinCopies)
                        .build())
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .minMinorAlleleCopies(somaticGeneCopyNumber.MinMinorAlleleCopyNumber)
                .build();
    }
}
