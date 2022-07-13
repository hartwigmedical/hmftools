package com.hartwig.hmftools.serve.sources.ckb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.ckb.classification.CkbClassificationConfig;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResourceTestFactory;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentApproachCurationType;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentApprochCurationEntry;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentAprroachCuration;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentAprroachCurationTest;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CkbExtractorTest {

    @Test
    public void canExtractFromCkbEntries() {
        EventClassifierConfig config = CkbClassificationConfig.build();
        CkbExtractor extractor = CkbExtractorFactory.buildCkbExtractor(config, RefGenomeResourceTestFactory.buildTestResource37());

        List<RelevantTreatmentApprochCurationEntry> curationEntries = Lists.newArrayList();
        curationEntries.add(RelevantTreatmentAprroachCurationTest.canGenerateCurationEntry(
                RelevantTreatmentApproachCurationType.TREATMENT_APPROACH_CURATION,
                "A",
                "A",
                "BRAF amplification",
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                "AA"));
        RelevantTreatmentAprroachCuration curator =
                new RelevantTreatmentAprroachCuration(curationEntries);

        List<CkbEntry> ckbEntries = Lists.newArrayList();
        ckbEntries.add(create("KIT", "amp", "KIT amp", "sensitive", "Actionable"));
        ckbEntries.add(create("BRAF", "V600E", "BRAF V600E", "sensitive", "Actionable"));
        ckbEntries.add(create("NTRK3", "fusion promiscuous", "NTRK3 fusion promiscuous", "sensitive", "Actionable"));
        ckbEntries.add(create("BRAF", "V600", "BRAF V600", "sensitive", "Actionable"));
        ckbEntries.add(create("BRAF", "exon 1 deletion", "BRAF exon 1 deletion", "sensitive", "Actionable"));
        ckbEntries.add(create("-", "MSI high", "MSI high", "sensitive", "Actionable"));
        ckbEntries.add(create("ALk", "EML4-ALK", "EML4-ALK Fusion", "sensitive", "Actionable"));

        ExtractionResult result = extractor.extract(ckbEntries, curator);
        assertEquals(1, result.knownHotspots().size());
        assertEquals(1, result.knownCopyNumbers().size());
        assertEquals(1, result.knownFusionPairs().size());
        assertEquals(1, result.actionableHotspots().size());
        assertEquals(2, result.actionableRanges().size());
        assertEquals(2, result.actionableGenes().size());
        assertEquals(1, result.actionableFusions().size());
        assertEquals(1, result.actionableCharacteristics().size());
    }

    @NotNull
    private static CkbEntry create(@NotNull String gene, @NotNull String variant, @NotNull String fullName, @NotNull String evidenceType,
            @NotNull String responseType) {
        return CkbTestFactory.createEntry(gene, variant, fullName, evidenceType, responseType, "AB", "cancer", "A", "DOID:162");
    }

    @Test
    public void canCurateCodons() {
        List<CodonAnnotation> codonAnnotations = Lists.newArrayList();
        CodonAnnotation codonAnnotation1 = ImmutableCodonAnnotation.builder()
                .gene("BRAF")
                .transcript("A")
                .chromosome("1")
                .start(10)
                .end(20)
                .mutationType(MutationTypeFilter.ANY)
                .rank(600)
                .build();

        CodonAnnotation codonAnnotation2 = ImmutableCodonAnnotation.builder()
                .gene("KRAS")
                .transcript("transcript")
                .chromosome("1")
                .start(10)
                .end(20)
                .mutationType(MutationTypeFilter.ANY)
                .rank(600)
                .build();

        codonAnnotations.add(codonAnnotation1);
        codonAnnotations.add(codonAnnotation2);

        List<CodonAnnotation> curatedCodons = CkbExtractor.curateCodons(codonAnnotations);

        CodonAnnotation codon1 = findByGene(curatedCodons, "BRAF");
        assertEquals(140753335, codon1.start());
        assertEquals(140753337, codon1.end());
        assertEquals("ENST00000646891", codon1.transcript());

        CodonAnnotation codon2 = findByGene(curatedCodons, "KRAS");
        assertEquals("KRAS", codon2.gene());
        assertEquals(10, codon2.start());
        assertEquals(20, codon2.end());
        assertEquals("transcript", codon2.transcript());
    }

    @NotNull
    private static CodonAnnotation findByGene(@NotNull Iterable<CodonAnnotation> codons, @NotNull String geneToFind) {
        for (CodonAnnotation codon : codons) {
            if (codon.gene().equals(geneToFind)) {
                return codon;
            }
        }

        throw new IllegalStateException("Could not find gene " + geneToFind);
    }
}