package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.algo.util.DriverUtils.convertReportedStatus;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.finding.ImmutableGainDeletion;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.orange.algo.util.DriverUtils;
import com.hartwig.hmftools.orange.report.finding.FindingKeys;

public class GermlineGainDeletionFactory
{
    private final EnsemblDataCache ensemblDataCache;

    public GermlineGainDeletionFactory(final EnsemblDataCache ensemblDataCache)
    {
        this.ensemblDataCache = ensemblDataCache;
    }

    public Map<GainDeletion, Boolean> getReportabilityMap(
            final List<GermlineDeletion> germlineDeletions, final List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        List<GermlineDeletion> germlineDeletionsHomozygousInTumor =
                germlineDeletions.stream().filter(d -> d.TumorStatus == GermlineStatus.HOM_DELETION).toList();
        Set<String> relevantGeneNames = germlineDeletionsHomozygousInTumor.stream().map(d -> d.GeneName).collect(Collectors.toSet());

        Map<GainDeletion, Boolean> delToReportability = Maps.newHashMap();
        for(String geneName : relevantGeneNames)
        {
            List<GermlineDeletion> deletionsForGene =
                    germlineDeletionsHomozygousInTumor.stream().filter(d -> d.GeneName.equals(geneName)).collect(Collectors.toList());
            GeneCopyNumber somaticGeneCopyNumber = GermlineDeletionUtil.findGeneCopyNumberForGene(geneName, allSomaticGeneCopyNumbers);

            GainDeletion del = toGainDel(geneName, deletionsForGene, somaticGeneCopyNumber);
            boolean reported = deletionsForGene.stream().anyMatch(d -> d.Reported == ReportedStatus.REPORTED);
            delToReportability.put(del, reported);
        }
        return delToReportability;
    }

    // FIX THIS: this is not correct
    private GainDeletion toGainDel(final String geneName, final List<GermlineDeletion> deletionsForGene,
            final GeneCopyNumber somaticGeneCopyNumber)
    {
        TranscriptData canonicalTranscript = GermlineDeletionUtil.findCanonicalTranscript(geneName, ensemblDataCache);
        CopyNumberInterpretation interpretation = GermlineDeletionUtil.deletionsCoverTranscript(deletionsForGene, canonicalTranscript)
                ? CopyNumberInterpretation.FULL_DEL
                : CopyNumberInterpretation.PARTIAL_DEL;
        double minCopies = GermlineDeletionUtil.getSomaticMinCopyNumber(deletionsForGene);
        double maxCopies = GermlineDeletionUtil.getSomaticMaxCopyNumber(deletionsForGene, somaticGeneCopyNumber, canonicalTranscript);
        String chromosome = GermlineDeletionUtil.getChromosome(deletionsForGene);
        String chromosomeBand = GermlineDeletionUtil.getChromosomeBand(deletionsForGene);

        // get the reported status
        com.hartwig.hmftools.datamodel.driver.ReportedStatus reportedStatus = DriverUtils.maxReportedStatus(
                deletionsForGene.stream().map(o -> convertReportedStatus(o.Reported)).collect(Collectors.toList())
        );

        // TODOHWL: double check this
        DriverInterpretation driverInterpretation = reportedStatus == com.hartwig.hmftools.datamodel.driver.ReportedStatus.REPORTED ?
                DriverInterpretation.HIGH : DriverInterpretation.LOW;

        return ImmutableGainDeletion.builder()
                .findingKey(FindingKeys.findingKey(geneName, interpretation, true))
                .reportedStatus(reportedStatus)
                .driverInterpretation(driverInterpretation)
                .interpretation(interpretation)
                .chromosome(chromosome)
                .chromosomeBand(chromosomeBand)
                .gene(geneName)
                .transcript(canonicalTranscript.TransName)
                .isCanonical(true)
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .build();
    }
}
