package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.Hotspot;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.PaveEntry;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class PurpleVariantFactory {

    private final PaveAlgo paveAlgo;

    @NotNull
    private final CompoundFilter mFilter;

    public PurpleVariantFactory withPassingOnlyFilter() {
        mFilter.add(new PassingVariantFilter());
        return this;
    }

    public PurpleVariantFactory(final PaveAlgo paveAlgo, final VariantContextFilter... filters) {
        this.paveAlgo = paveAlgo;
        mFilter = new CompoundFilter(true);
        mFilter.addAll(Arrays.asList(filters));
        mFilter.add(new HumanChromosomeFilter());
        mFilter.add(new NTFilter());
    }

    public List<PurpleVariant> fromVCFFile(final String tumor, @Nullable final String reference, @Nullable final String rna,
            final String vcfFile) throws IOException {
        List<PurpleVariant> result = new ArrayList<>();

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false)) {
            final VCFHeader header = (VCFHeader) reader.getHeader();

            if (!sampleInFile(tumor, header)) {
                throw new IllegalArgumentException("Sample " + tumor + " not found in vcf file " + vcfFile);
            }

            if (reference != null && !sampleInFile(reference, header)) {
                throw new IllegalArgumentException("Sample " + reference + " not found in vcf file " + vcfFile);
            }

            if (rna != null && !sampleInFile(rna, header)) {
                throw new IllegalArgumentException("Sample " + rna + " not found in vcf file " + vcfFile);
            }

            if (!header.hasFormatLine("AD")) {
                throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
            }

            for (VariantContext variant : reader.iterator()) {
                if (mFilter.test(variant)) {
                    try {
                        PurpleVariant purpleSomaticVariant = createVariant(tumor, reference, rna, variant);
                        result.add(purpleSomaticVariant);
                    } catch (IllegalArgumentException e) {
                        // ignore, consider the sample filtered
                    }
                }
            }
        }

        return result;
    }

    public PurpleVariant createVariant(final String sample, @Nullable final String reference, @Nullable final String rna,
            final VariantContext context) {
        if (!mFilter.test(context)) {
            throw new IllegalArgumentException(String.format("Variant could not be created because sample [%s] does not have status PASS", sample));
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

        return helperCreateVariant(context, tumorDepth, reference, rna);
    }

    private PurpleVariant helperCreateVariant(VariantContext context, AllelicDepth tumorDepth, @Nullable String reference,
            @Nullable String rna) {
        VariantContextDecorator contextDecorator = new VariantContextDecorator(context);
        final VariantImpact variantImpact = contextDecorator.variantImpact();
        final List<VariantTranscriptImpact> variantTranscriptImpacts = VariantTranscriptImpact.fromVariantContext(context);
        final List<PurpleVariantTranscriptImpact> purpleVariantTranscriptImpacts =
                variantTranscriptImpacts.stream().map(PurpleConversion::convert).collect(Collectors.toList());
        final Optional<AllelicDepth> rnaDepth = rnaDepth(context, rna);
        final GenotypeStatus genotypeStatus = reference != null ? contextDecorator.genotypeStatus(reference) : null;

        return ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.valueOf(contextDecorator.type().name()))
                .gene(variantImpact.CanonicalGeneName)
                .chromosome(contextDecorator.chromosome())
                .position(contextDecorator.position())
                .ref(contextDecorator.ref())
                .alt(contextDecorator.alt())
                .worstCodingEffect(PurpleConversion.convert(variantImpact.WorstCodingEffect))
                .canonicalImpact(extractCanonicalImpact(contextDecorator))
                .otherImpacts(extractOtherImpacts(contextDecorator))
                .hotspot(Hotspot.valueOf(contextDecorator.hotspot().name()))
                .reported(contextDecorator.reported())
                .tumorDepth(extractTumorDepth(tumorDepth))
                .rnaDepth(rnaDepth.map(PurpleConversion::convert).orElse(null))
                .adjustedCopyNumber(contextDecorator.adjustedCopyNumber())
                .adjustedVAF(contextDecorator.adjustedVaf())
                .minorAlleleCopyNumber(contextDecorator.minorAlleleCopyNumber())
                .variantCopyNumber(contextDecorator.variantCopyNumber())
                .biallelic(contextDecorator.biallelic())
                .genotypeStatus(PurpleGenotypeStatus.valueOf((genotypeStatus != null ? genotypeStatus : GenotypeStatus.UNKNOWN).name()))
                .repeatCount(contextDecorator.repeatCount())
                .subclonalLikelihood(context.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0))
                .localPhaseSets(context.getAttributeAsIntList(LOCAL_PHASE_SET, 0))
                .variantTranscriptImpacts(purpleVariantTranscriptImpacts)
                .build();
    }

    private PurpleTranscriptImpact extractCanonicalImpact(VariantContextDecorator contextDecorator) {
        var variantImpact = contextDecorator.variantImpact();

        PaveEntry paveEntry = paveAlgo.run(variantImpact.CanonicalGeneName, variantImpact.CanonicalTranscript, contextDecorator.position());

        List<VariantEffect> variantEffects = VariantEffect.effectsToList(variantImpact.CanonicalEffect);
        List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapToList(variantEffects, PurpleConversion::convert);
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(variantImpact.CanonicalTranscript)
                .hgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .hgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .effects(purpleVariantEffects)
                .codingEffect(PurpleConversion.convert(variantImpact.CanonicalCodingEffect))
                .build();
    }

    private List<PurpleTranscriptImpact> extractOtherImpacts(VariantContextDecorator contextDecorator) {

        var variantImpact = contextDecorator.variantImpact();
        PaveEntry paveEntry = paveAlgo.run(variantImpact.CanonicalGeneName, variantImpact.CanonicalTranscript, contextDecorator.position());

        List<PurpleTranscriptImpact> otherImpacts = Lists.newArrayList();
        for (AltTranscriptReportableInfo altInfo : AltTranscriptReportableInfo.parseAltTranscriptInfo(variantImpact.OtherReportableEffects)) {
            List<VariantEffect> variantEffects = VariantEffect.effectsToList(altInfo.Effects);
            List<PurpleVariantEffect> purpleVariantEffects = ConversionUtil.mapToList(variantEffects, PurpleConversion::convert);
            otherImpacts.add(ImmutablePurpleTranscriptImpact.builder()
                    .transcript(altInfo.TransName)
                    .hgvsCodingImpact(altInfo.HgvsCoding)
                    .hgvsProteinImpact(altInfo.HgvsProtein)
                    .affectedCodon(paveEntry != null ? paveEntry.affectedCodon() : null)
                    .affectedExon(paveEntry != null ? paveEntry.affectedExon() : null)
                    .spliceRegion(variantImpact.CanonicalSpliceRegion)
                    .effects(purpleVariantEffects)
                    .codingEffect(PurpleConversion.convert(altInfo.Effect))
                    .build());
        }
        return otherImpacts;
    }

    private static PurpleAllelicDepth extractTumorDepth(AllelicDepth tumorDepth) {
        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(tumorDepth.alleleReadCount())
                .totalReadCount(tumorDepth.totalReadCount())
                .build();
    }

    private static Optional<AllelicDepth> rnaDepth(VariantContext context, String rna) {
        return Optional.ofNullable(context.getGenotype(rna)).filter(AllelicDepth::containsAllelicDepth).map(AllelicDepth::fromGenotype);
    }

    private static boolean sampleInFile(final String sample, final VCFHeader header) {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }
}