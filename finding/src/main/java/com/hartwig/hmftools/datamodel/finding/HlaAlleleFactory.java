package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

// lilac shows two copies of HLA alleles even if they are the
// same allele. This class combine the alleles and sum the copies
public class HlaAlleleFactory
{
    private static final Logger LOGGER = LogManager.getLogger(HlaAlleleFactory.class);

    private static final String PASS = "PASS";

    private HlaAlleleFactory() {
    }

    public static Findings<HlaAllele> createHlaAllelesFindings(@NotNull OrangeRecord orangeRecord, boolean hasReliablePurity)
    {
        LilacRecord lilac = orangeRecord.lilac();
        return ImmutableFindings.<HlaAllele>builder()
                .status(lilac != null ? lilac.qc().equals(PASS) ? FindingsStatus.OK : FindingsStatus.NOT_RELIABLE : FindingsStatus.NOT_AVAILABLE)
                .all(lilac != null ? HlaAlleleFactory.convertHlaAlleles(lilac, hasReliablePurity, !orangeRecord.tumorOnlyMode(),
                        orangeRecord.isofox() != null) : List.of())
                .build();
    }

    @NotNull
    public static List<HlaAllele> convertHlaAlleles(@NotNull LilacRecord lilac, boolean hasReliablePurity, boolean hasRef, boolean hasRna)
    {
        Map<String, List<LilacAllele>> hlaAllelesMap = lilac.alleles()
                .stream()
                .collect(Collectors.groupingBy(LilacAllele::allele));

        List<HlaAllele> hlaAlleles = new ArrayList<>();
        for (Map.Entry<String, List<LilacAllele>> keyMap : hlaAllelesMap.entrySet())
        {
            LilacAllele lilacAllele = keyMap.getValue().get(0);

            // NOTE: the fragment counts are doubled in lilac if an allele is present twice
            ImmutableHlaAllele.Builder builder = ImmutableHlaAllele.builder()
                    .findingKey(FindingKeys.hlaAllele(lilacAllele))
                    .allele(lilacAllele.allele())
                    .gene(extractHLAGene(lilacAllele.allele()))
                    .refFragments(hasRef ? lilacAllele.refFragments() : null)
                    .tumorFragments(lilacAllele.tumorFragments())
                    .rnaFragments(hasRna ? lilacAllele.rnaFragments() : null)
                    .somaticMissense(lilacAllele.somaticMissense())
                    .somaticNonsenseOrFrameshift(lilacAllele.somaticNonsenseOrFrameshift())
                    .somaticSplice(lilacAllele.somaticSplice())
                    .somaticSynonymous(lilacAllele.somaticSynonymous())
                    .somaticInframeIndel(lilacAllele.somaticInframeIndel());

            if (keyMap.getValue().size() == 1)
            {
                hlaAlleles.add(builder
                        .germlineCopyNumber(1)
                        .tumorCopyNumber(hasReliablePurity ? lilacAllele.tumorCopyNumber() : null)
                        .build());

            } else if (keyMap.getValue().size() == 2)
            {
                LilacAllele allele2 = keyMap.getValue().get(1);
                double tumorCopies = lilacAllele.tumorCopyNumber() + allele2.tumorCopyNumber();

                hlaAlleles.add(builder
                        .germlineCopyNumber(2)
                        .tumorCopyNumber(hasReliablePurity ? tumorCopies : null)
                        .build());
            } else
            {
                LOGGER.warn("To many hla alleles of allele '{}'", keyMap.getKey());
            }
        }

        hlaAlleles.sort(HlaAllele.COMPARATOR);
        return hlaAlleles;
    }

    @NotNull
    public static String extractHLAGene(@NotNull String allele) {
        if (allele.startsWith("A*")) {
            return "HLA-A";
        } else if (allele.startsWith("B*")) {
            return "HLA-B";
        } else if (allele.startsWith("C*")) {
            return "HLA-C";
        } else {
            LOGGER.warn("Unknown HLA gene name '{}' present! ", allele);
            return Strings.EMPTY;
        }
    }
}
