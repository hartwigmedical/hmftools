package com.hartwig.hmftools.patientdb.readers;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.consensus.ConsensusRule;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientdb.data.SomaticVariantData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SomaticVariantReader {
    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantReader.class);

    private static final String SOMATIC_EXTENSION = "_melted.vcf";

    @NotNull
    private final ConsensusRule consensusRule;

    private static final List<VariantConsequence> ACTIONABLE_CONSEQUENCES = Lists.newArrayList(
            VariantConsequence.TRANSCRIPT_ABLATION, VariantConsequence.TRANSCRIPT_AMPLIFICATION,
            VariantConsequence.SPLICE_ACCEPTOR_VARIANT, VariantConsequence.SPLICE_DONOR_VARIANT,
            VariantConsequence.SPLICE_REGION_VARIANT, VariantConsequence.STOP_GAINED, VariantConsequence.STOP_LOST,
            VariantConsequence.START_LOST, VariantConsequence.FRAMESHIFT_VARIANT, VariantConsequence.INFRAME_INSERTION,
            VariantConsequence.INFRAME_DELETION, VariantConsequence.MISSENSE_VARIANT);

    public SomaticVariantReader(@NotNull final ConsensusRule consensusRule) throws IOException, HartwigException {
        this.consensusRule = consensusRule;
    }

    @NotNull
    public List<SomaticVariantData> read(@NotNull final String runDirectory) throws IOException, HartwigException {
        final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(runDirectory, SOMATIC_EXTENSION);
        final List<SomaticVariant> passedVariants = passOnly(variantFile.variants());
        final List<SomaticVariant> consensusPassedVariants = consensusRule.removeUnreliableVariants(passedVariants);
        final List<SomaticVariant> actionableVariants = Lists.newArrayList();
        consensusPassedVariants.forEach(consensusPassedVariant -> consensusPassedVariant.annotations().forEach(
                annotation -> annotation.consequences().forEach(consequence -> {
                    if (ACTIONABLE_CONSEQUENCES.contains(consequence)) {
                        if (!actionableVariants.contains(consensusPassedVariant)) {
                            actionableVariants.add(consensusPassedVariant);
                        }
                    }
                })));
        return toVariantData(actionableVariants);
    }

    @NotNull
    private static List<SomaticVariantData> toVariantData(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariantData> somaticVariantDatas = Lists.newArrayList();
        for (final SomaticVariant variant : variants) {
            // MIVO: get the first annotation for now, eventually we will want all
            final VariantAnnotation variantAnnotation = variant.annotations().get(0);
            variant.annotations().forEach(annotation -> {
                if (!annotation.gene().equals(variantAnnotation.gene())) {
                    LOGGER.warn("Annotated gene (" + annotation.gene()
                            + ") does not match gene expected from first annotation ( " + variantAnnotation.gene()
                            + ") for variant: " + variant);
                }
            });
            final SomaticVariantData somaticVariantData = new SomaticVariantData(variantAnnotation.gene(),
                    variant.chromosome() + ":" + variant.position(), variant.ref(), variant.alt(), variant.cosmicID(),
                    variant.totalReadCount(), variant.alleleReadCount());
            somaticVariantDatas.add(somaticVariantData);
        }
        return somaticVariantDatas;
    }
}
