package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

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

public class PurpleVariantContextLoader
{
    @NotNull
    private final CompoundFilter filter;

    public static PurpleVariantContextLoader withPassingOnlyFilter()
    {
        return new PurpleVariantContextLoader(new PassingVariantFilter());
    }

    public PurpleVariantContextLoader(final VariantContextFilter... filters)
    {
        filter = new CompoundFilter(true);
        filter.addAll(Arrays.asList(filters));
        filter.add(new HumanChromosomeFilter());
        filter.add(new NTFilter());
    }

    public List<PurpleVariantContext> fromVCFFile(final String tumor, @Nullable final String reference, @Nullable final String rna,
            final String vcfFile) throws IOException
    {
        List<PurpleVariantContext> result = new ArrayList<>();

        final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);
        final VCFHeader header = (VCFHeader) reader.getHeader();

        if(!sampleInFile(tumor, header))
        {
            throw new IllegalArgumentException("Sample " + tumor + " not found in vcf file " + vcfFile);
        }

        if(reference != null && !sampleInFile(reference, header))
        {
            throw new IllegalArgumentException("Sample " + reference + " not found in vcf file " + vcfFile);
        }

        if(rna != null && !sampleInFile(rna, header))
        {
            throw new IllegalArgumentException("Sample " + rna + " not found in vcf file " + vcfFile);
        }

        if(!header.hasFormatLine("AD"))
        {
            throw new IllegalArgumentException("Allelic depths is a required format field in vcf file " + vcfFile);
        }

        for(VariantContext variantContext : reader.iterator())
        {
            if(filter.test(variantContext))
            {
                PurpleVariantContext purpleVariantContext = createPurpleVariantContext(variantContext, tumor, reference, rna);
                result.add(purpleVariantContext);
            }
        }
        return result;
    }

    public PurpleVariantContext createPurpleVariantContext(VariantContext variantContext, String sample, @Nullable String reference,
            @Nullable String rna)
    {
        if(!AllelicDepth.containsAllelicDepth(variantContext.getGenotype(sample)))
        {
            throw new IllegalArgumentException(String.format(
                    "Variant could not be created because sample [%s] does not contain allelic depth",
                    sample));
        }

        final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(variantContext.getGenotype(sample));
        int readCount = tumorDepth.totalReadCount();

        if(readCount <= 0)
        {
            throw new IllegalArgumentException(String.format(
                    "Variant could not be created because tumor depth read count should be greater than 0 (actual value: %s)",
                    readCount));
        }

        return createCreatePurpleVariantContext(variantContext, tumorDepth, reference, rna);
    }

    private PurpleVariantContext createCreatePurpleVariantContext(VariantContext variantContext, AllelicDepth tumorDepth,
            @Nullable String reference, @Nullable String rna)
    {
        VariantContextDecorator contextDecorator = new VariantContextDecorator(variantContext);
        final VariantImpact variantImpact = contextDecorator.variantImpact();
        final List<VariantTranscriptImpact> variantTranscriptImpacts = VariantTranscriptImpact.fromVariantContext(variantContext);
        final AllelicDepth rnaDepth = extractRnaDepth(variantContext, rna);

        return ImmutablePurpleVariantContext.builder()
                .chromosome(contextDecorator.chromosome())
                .position(contextDecorator.position())
                .totalReadCount(tumorDepth.totalReadCount())
                .alleleReadCount(tumorDepth.alleleReadCount())
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .type(contextDecorator.type())
                .gene(variantImpact.CanonicalGeneName)
                .ref(contextDecorator.ref())
                .alt(contextDecorator.alt())
                .canonicalTranscript(variantImpact.CanonicalTranscript)
                .canonicalEffect(variantImpact.CanonicalEffect)
                .canonicalCodingEffect(variantImpact.CanonicalCodingEffect)
                .canonicalHgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .canonicalHgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .worstCodingEffect(variantImpact.WorstCodingEffect)
                .otherImpacts(variantTranscriptImpacts)
                .hotspot(contextDecorator.hotspot())
                .reported(contextDecorator.reported())
                .tumorDepth(tumorDepth)
                .rnaDepth(rnaDepth)
                .adjustedCopyNumber(contextDecorator.adjustedCopyNumber())
                .adjustedVAF(contextDecorator.adjustedVaf())
                .minorAlleleCopyNumber(contextDecorator.minorAlleleCopyNumber())
                .variantCopyNumber(contextDecorator.variantCopyNumber())
                .biallelic(contextDecorator.biallelic())
                .genotypeStatus(contextDecorator.genotypeStatus(reference))
                .repeatCount(contextDecorator.repeatCount())
                .subclonalLikelihood(variantContext.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0))
                .localPhaseSets(variantContext.getAttributeAsIntList(LOCAL_PHASE_SET, 0))
                .build();
    }

    @Nullable
    private static AllelicDepth extractRnaDepth(VariantContext context, @Nullable String rna)
    {
        return Optional.ofNullable(context.getGenotype(rna))
                .filter(AllelicDepth::containsAllelicDepth)
                .map(AllelicDepth::fromGenotype)
                .orElse(null);
    }

    private static boolean sampleInFile(final String sample, final VCFHeader header)
    {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }
}
