package com.hartwig.hmftools.finding;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.HlaAlleleBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// lilac shows two copies of HLA alleles even if they are the
// same allele. This class combine the alleles and sum the copies
public class HlaAlleleFactory
{
    private static final Logger LOGGER = LogManager.getLogger(HlaAlleleFactory.class);

    private static final Pattern HLA_REGEX = Pattern.compile("""
            ^(?<geneSymbol>\\w+)\\*(?<alleleGroup>\\d{2}):(?<hlaProtein>\\d{2,3})N?$""");
    private static final String PASS = "PASS";

    private HlaAlleleFactory()
    {
    }

    public static FindingList<HlaAllele> createHlaAllelesFindings(OrangeRecord orangeRecord, FindingStatus findingStatus)
    {
        LilacRecord lilac = orangeRecord.lilac();
        if(lilac != null)
        {
            return FindingListBuilder.<HlaAllele>builder()
                    .status(findingStatus)
                    .findings(HlaAlleleFactory.convertHlaAlleles(lilac,
                            !orangeRecord.tumorOnlyMode(),
                            orangeRecord.isofox() != null))
                    .build();
        }
        else
        {
            throw new IllegalStateException("Missing Lilac record in Orange record");
        }
    }

    public static List<HlaAllele> convertHlaAlleles(LilacRecord lilac, boolean hasRef, boolean hasRna)
    {
        Map<String, Integer> hlaAlleleCount = new HashMap<>();

        return lilac.alleles().stream()
                .map(lilacAllele -> convertLilacAllele(lilacAllele, hasRef, hasRna,
                        hlaAlleleCount.compute(lilacAllele.allele(),
                                (key, oldValue) -> oldValue == null ? 1 : oldValue + 1)))
                .sorted(HlaAllele.COMPARATOR)
                .toList();
    }

    static HlaAllele convertLilacAllele(LilacAllele lilacAllele, boolean hasRef, boolean hasRna, int alleleCopy)
    {
        var matcher = matchHlaRegEx(lilacAllele.allele());
        String gene = matcher.group("geneSymbol");
        HlaAllele.GeneClass geneClass = lilacAllele.geneClass().equals(HlaCommon.MHC_CLASS_I) ? HlaAllele.GeneClass.MHC_CLASS_I : HlaAllele.GeneClass.MHC_CLASS_II;
        String alleleGroup = matcher.group("alleleGroup");
        String hlaProtein = matcher.group("hlaProtein");

        HlaAlleleBuilder builder = HlaAlleleBuilder.builder()
                .findingKey(FindingKeys.hlaAllele(lilacAllele, alleleCopy))
                .geneClass(geneClass)
                .geneSymbol(gene)
                .alleleGroup(alleleGroup)
                .hlaProtein(hlaProtein)
                .qcStatus(Arrays.stream(lilacAllele.qcStatus().split(";"))
                        .map(HlaAllele.QcStatus::valueOf)
                        .collect(Collectors.toSet()))
                .germlineCopyNumber(1)
                .tumorCopyNumber(lilacAllele.tumorCopyNumber())
                .refFragments(hasRef ? lilacAllele.refFragments() : null)
                .tumorFragments(lilacAllele.tumorFragments())
                .rnaFragments(hasRna ? lilacAllele.rnaFragments() : null)
                .somaticMissense(lilacAllele.somaticMissense())
                .somaticNonsenseOrFrameshift(lilacAllele.somaticNonsenseOrFrameshift())
                .somaticSplice(lilacAllele.somaticSplice())
                .somaticSynonymous(lilacAllele.somaticSynonymous())
                .somaticInframeIndel(lilacAllele.somaticInframeIndel());

        return builder.build();
    }

    @NotNull
    static Matcher matchHlaRegEx(String allele)
    {
        var matcher = HLA_REGEX.matcher(allele);
        if(!matcher.matches())
        {
            throw new IllegalStateException("Can't extract HLA gene, alleleGroup and hlaProtein from " + allele);
        }
        return matcher;
    }
}
