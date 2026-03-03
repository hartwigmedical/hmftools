package com.hartwig.hmftools.finding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.hartwig.hmftools.finding.datamodel.FindingList;
import com.hartwig.hmftools.finding.datamodel.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.HlaAlleleBuilder;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// lilac shows two copies of HLA alleles even if they are the
// same allele. This class combine the alleles and sum the copies
public class HlaAlleleFactory
{
    private static final Logger LOGGER = LogManager.getLogger(HlaAlleleFactory.class);

    private static final Pattern HLA_REGEX = Pattern.compile("""
            ^(?<gene>[A-Z]+)\\*(?<alleleGroup>\\d{2}):(?<hlaProtein>\\d{2,3})N?$""");
    private static final String PASS = "PASS";

    private HlaAlleleFactory()
    {
    }

    public static FindingList<HlaAllele> createHlaAllelesFindings(OrangeRecord orangeRecord, boolean hasReliablePurity,
            EventFactory eventFactory)
    {
        LilacRecord lilac = orangeRecord.lilac();
        if(lilac != null)
        {
            return FindingListBuilder.<HlaAllele>builder()
                    .status(lilac.qc().equals(PASS) ? FindingsStatus.OK : FindingsStatus.NOT_RELIABLE)
                    .findings(HlaAlleleFactory.convertHlaAlleles(lilac,
                            hasReliablePurity,
                            !orangeRecord.tumorOnlyMode(),
                            orangeRecord.isofox() != null, eventFactory))
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingList();
        }
    }

    public static List<HlaAllele> convertHlaAlleles(LilacRecord lilac, boolean hasReliablePurity, boolean hasRef, boolean hasRna,
            EventFactory eventFactory)
    {
        Map<String, List<LilacAllele>> hlaAllelesMap = lilac.alleles()
                .stream()
                .collect(Collectors.groupingBy(LilacAllele::allele));

        // newer version of Lilac puts acStatus in alleles instead. This backport version will do the same
        Set<HlaAllele.QcStatus> qcStatus = Arrays.stream(lilac.qc().split(";"))
                .map(HlaAllele.QcStatus::valueOf)
                .collect(Collectors.toSet());

        List<HlaAllele> hlaAlleles = new ArrayList<>();
        for(Map.Entry<String, List<LilacAllele>> keyMap : hlaAllelesMap.entrySet())
        {
            LilacAllele lilacAllele = keyMap.getValue().get(0);

            var matcher = HLA_REGEX.matcher(lilacAllele.allele());
            //throw IllegalStateException("Can't extract HLA gene, alleleGroup and hlaProtein from ${allele.allele()}")
            String gene = matcher.group("gene");
            String geneClass = "HLA_" + gene;
            String alleleGroup = matcher.group("alleleGroup");
            String hlaProtein = matcher.group("hlaProtein");

            // NOTE: the fragment counts are doubled in lilac if an allele is present twice
            HlaAlleleBuilder builder = HlaAlleleBuilder.builder()
                    .findingKey(FindingKeys.hlaAllele(lilacAllele))
                    .event(eventFactory.immunologyEvent(lilacAllele))
                    .geneClass(geneClass)
                    .gene(gene)
                    .allele(lilacAllele.allele())
                    .alleleGroup(alleleGroup)
                    .hlaProtein(hlaProtein)
                    .qcStatus(qcStatus)
                    .refFragments(hasRef ? lilacAllele.refFragments() : null)
                    .tumorFragments(lilacAllele.tumorFragments())
                    .rnaFragments(hasRna ? lilacAllele.rnaFragments() : null)
                    .somaticMissense(lilacAllele.somaticMissense())
                    .somaticNonsenseOrFrameshift(lilacAllele.somaticNonsenseOrFrameshift())
                    .somaticSplice(lilacAllele.somaticSplice())
                    .somaticSynonymous(lilacAllele.somaticSynonymous())
                    .somaticInframeIndel(lilacAllele.somaticInframeIndel());

            if(keyMap.getValue().size() == 1)
            {
                hlaAlleles.add(builder
                        .germlineCopyNumber(1)
                        .tumorCopyNumber(hasReliablePurity ? lilacAllele.tumorCopyNumber() : null)
                        .build());

            }
            else if(keyMap.getValue().size() == 2)
            {
                LilacAllele allele2 = keyMap.getValue().get(1);
                double tumorCopies = lilacAllele.tumorCopyNumber() + allele2.tumorCopyNumber();

                hlaAlleles.add(builder
                        .germlineCopyNumber(2)
                        .tumorCopyNumber(hasReliablePurity ? tumorCopies : null)
                        .build());
            }
            else
            {
                LOGGER.warn("To many hla alleles of allele '{}'", keyMap.getKey());
            }
        }

        hlaAlleles.sort(HlaAllele.COMPARATOR);
        return hlaAlleles;
    }
}
