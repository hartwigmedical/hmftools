package com.hartwig.hmftools.orange.algo.immuno;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.orange.OrangeConfig;

public final class LilacInterpreter
{
    public static LilacRecord build(final OrangeConfig config) throws IOException
    {
        String qcFile = LilacQcData.generateFilename(config.LilacDir, config.TumorId);
        List<LilacQcData> qcData = LilacQcData.read(qcFile);

        String allelesFile = LilacAllele.generateFilename(config.LilacDir, config.TumorId);
        List<LilacAllele> alleles = LilacAllele.read(allelesFile);

        List<com.hartwig.hmftools.datamodel.hla.LilacAllele> convertedAlleles = Lists.newArrayList();

        for(LilacAllele allele : alleles)
        {
            LilacQcData alleleQcData = qcData.stream().filter(x -> x.genes().equals(allele.genes())).findFirst().orElse(null);
            convertedAlleles.add(convert(allele, alleleQcData.status(), config.hasReference(), config.hasRNA()));
        }

        return ImmutableLilacRecord.builder().alleles(convertedAlleles).build();
    }

    private static com.hartwig.hmftools.datamodel.hla.LilacAllele convert(
            final com.hartwig.hmftools.common.hla.LilacAllele allele, final String qcStatus,  boolean hasRef, boolean hasRna)
    {
        return ImmutableLilacAllele.builder()
                .geneClass(allele.genes())
                .allele(allele.allele())
                .qcStatus(qcStatus)
                .tumorCopyNumber(allele.tumorCopyNumber())
                .refFragments(hasRef ? allele.refFragments() : null)
                .tumorFragments(allele.tumorFragments())
                .rnaFragments(hasRna ? allele.rnaFragments() : null)
                .somaticMissense(allele.somaticMissense())
                .somaticNonsenseOrFrameshift(allele.somaticNonsenseOrFrameshift())
                .somaticSplice(allele.somaticSplice())
                .somaticSynonymous(allele.somaticSynonymous())
                .somaticInframeIndel(allele.somaticInframeIndel())
                .build();
    }

}
