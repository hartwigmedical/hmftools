package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class PurpleSomaticVariantFactory {

    private final PassingVariantFilter passingVariantFilter = new PassingVariantFilter();
    private int mCreatedCount;
    private int mFilteredCount;

    private static final String RECOVERED_FLAG = "RECOVERED";

    public int getCreatedCount() {
        return mCreatedCount;
    }

    public int getFilteredCount() {
        return mFilteredCount;
    }

    public List<PurpleSomaticVariant> fromVCFFile(final String tumor, @Nullable final String reference, @Nullable final String rna,
            final String vcfFile, boolean useCheckReference) throws IOException {
        List<PurpleSomaticVariant> result = new ArrayList<>();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {
            final VCFHeader header = (VCFHeader) reader.getHeader();

            if (!sampleInFile(tumor, header)) {
                throw new IllegalArgumentException("Sample " + tumor + " not found in vcf file " + vcfFile);
            }

            if (useCheckReference && reference != null && !sampleInFile(reference, header)) {
                throw new IllegalArgumentException("Sample " + reference + " not found in vcf file " + vcfFile);
            }

            if (rna != null && !sampleInFile(rna, header)) {
                throw new IllegalArgumentException("Sample " + rna + " not found in vcf file " + vcfFile);
            }

            if (!header.hasFormatLine("AD")) {
                throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
            }

            for (VariantContext variant : reader.iterator()) {
                if (passingVariantFilter.test(variant)) {
                    try {
                        PurpleSomaticVariant purpleSomaticVariant = createVariant(tumor, reference, rna, variant);
                        result.add(purpleSomaticVariant);
                        ++mCreatedCount;
                    } catch (IllegalArgumentException e) {
                        ++mFilteredCount;
                    }
                } else {
                    ++mFilteredCount;
                }
            }
        }

        return result;
    }

    //TODO find out if exceptions are appropriate or if optional is preferred
    public PurpleSomaticVariant createVariant(final String sample, @Nullable final String reference, @Nullable final String rna,
            final VariantContext context) {
        if (!passingVariantFilter.test(context)) {
            throw new IllegalArgumentException(String.format("Variant could not be created because sample [%s] does not PASS", sample));
        }

        if (!AllelicDepth.containsAllelicDepth(context.getGenotype(sample))) {
            throw new IllegalArgumentException(String.format(
                    "Variant could not be created because sample [%s] does not contain allelic depth",
                    sample));
        }

        final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));
        int readCount = tumorDepth.totalReadCount();

        if (readCount <= 0) {
            throw new IllegalArgumentException(String.format(
                    "Variant could not be created because tumor depth read count should be greater than 0 (actual value: %s)",
                    readCount));
        }

        return helperCreateVariant(reference, rna, tumorDepth, context);
    }

    private PurpleSomaticVariant helperCreateVariant(String reference, String rna, AllelicDepth tumorDepth, VariantContext context) {
        VariantContextDecorator contextDecorator = new VariantContextDecorator(context);
        final VariantImpact variantImpact = contextDecorator.variantImpact();
        final List<VariantTranscriptImpact> variantTranscriptImpacts = VariantTranscriptImpact.fromVariantContext(context);

        final Optional<AllelicDepth> referenceDepth = Optional.ofNullable(reference)
                .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                .filter(AllelicDepth::containsAllelicDepth)
                .map(AllelicDepth::fromGenotype);

        final Optional<AllelicDepth> rnaDepth = Optional.ofNullable(rna)
                .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                .filter(AllelicDepth::containsAllelicDepth)
                .map(AllelicDepth::fromGenotype);

        return ImmutablePurpleSomaticVariantImpl.builder()
                .qual(contextDecorator.qual())
                .type(contextDecorator.type())
                .filter(contextDecorator.filter())
                .chromosome(contextDecorator.chromosome())
                .position(contextDecorator.position())
                .ref(contextDecorator.ref())
                .alt(contextDecorator.alt())
                .alleleReadCount(tumorDepth.alleleReadCount())
                .totalReadCount(tumorDepth.totalReadCount())
                .hotspot(contextDecorator.hotspot())
                .minorAlleleCopyNumber(contextDecorator.minorAlleleCopyNumber())
                .adjustedCopyNumber(contextDecorator.adjustedCopyNumber())
                .adjustedVAF(contextDecorator.adjustedVaf())
                .variantCopyNumber(contextDecorator.variantCopyNumber())
                .mappability(contextDecorator.mappability())
                .tier(contextDecorator.tier())
                .trinucleotideContext(contextDecorator.trinucleotideContext())
                .microhomology(contextDecorator.microhomology())
                .repeatCount(contextDecorator.repeatCount())
                .repeatSequence(contextDecorator.repeatSequence())
                .reported(contextDecorator.reported())
                .biallelic(contextDecorator.biallelic())
                .gene(variantImpact.CanonicalGeneName)
                .canonicalTranscript(variantImpact.CanonicalTranscript)
                .canonicalEffect(variantImpact.CanonicalEffect)
                .canonicalCodingEffect(variantImpact.CanonicalCodingEffect)
                .canonicalHgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .canonicalHgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .otherReportedEffects(variantImpact.OtherReportableEffects)
                .worstCodingEffect(variantImpact.WorstCodingEffect)
                .genesAffected(variantImpact.GenesAffected)
                .subclonalLikelihood(context.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0))
                .germlineStatus(GermlineStatus.valueOf(context.getAttributeAsString(PURPLE_GERMLINE_INFO, "UNKNOWN")))
                .kataegis(context.getAttributeAsString(KATAEGIS_FLAG, Strings.EMPTY))
                .recovered(context.getAttributeAsBoolean(RECOVERED_FLAG, false))
                .clinvarInfo(context.getAttributeAsString(CLNSIG, ""))
                .gnomadFrequency(context.getAttributeAsDouble(GNOMAD_FREQ, 0))
                .somaticLikelihood(SomaticLikelihood.valueOf(context.getAttributeAsString(PANEL_SOMATIC_LIKELIHOOD,
                        SomaticLikelihood.UNKNOWN.toString())))
                .variantTranscriptImpacts(variantTranscriptImpacts)
                .localPhaseSets(context.getAttributeAsIntList(LOCAL_PHASE_SET, 0))
                .referenceDepth(referenceDepth.orElse(null))
                .rnaDepth(rnaDepth.orElse(null))
                .build();
    }

    private static boolean sampleInFile(final String sample, final VCFHeader header) {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

}