package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFHeader;

public class SmallVariantFactory implements VariantContextFilter
{
    public static SmallVariantFactory passOnlyInstance()
    {
        return new SmallVariantFactory(new PassingVariantFilter());
    }

    private final CompoundFilter mFilter;
    private int mCreatedCount;
    private int mFilteredCount;

    private boolean mDropDuplicates;
    private SmallVariant mLastNonFilteredVariant;

    public SmallVariantFactory(final VariantContextFilter... filters)
    {
        mFilter = new CompoundFilter(true);
        mFilter.addAll(Arrays.asList(filters));
        mFilter.add(new HumanChromosomeFilter());
        mFilter.add(new NTFilter());
        mCreatedCount = 0;
        mFilteredCount = 0;
        mLastNonFilteredVariant = null;
        mDropDuplicates = false;
    }

    public int getCreatedCount()
    {
        return mCreatedCount;
    }

    public int getFilteredCount()
    {
        return mFilteredCount;
    }

    public void setDropDuplicates() { mDropDuplicates = true; }

    public static List<SmallVariant> loadVariants(final String tumor, final String vcfFile) throws IOException
    {
        return new SmallVariantFactory().fromVCFFile(tumor, null, null, vcfFile);
    }

    public List<SmallVariant> fromVCFFile(final String tumor, final String vcfFile) throws IOException
    {
        return fromVCFFile(tumor, null, null, vcfFile);
    }

    public List<SmallVariant> fromVCFFile(final String tumor, @Nullable final String reference,
            final String vcfFile) throws IOException
    {
        final List<SmallVariant> result = Lists.newArrayList();
        fromVCFFile(tumor, reference, null, vcfFile, true, result::add);
        return result;
    }

    public List<SmallVariant> fromVCFFile(
            final String tumor, @Nullable final String reference, @Nullable final String rna, final String vcfFile) throws IOException
    {
        final List<SmallVariant> result = Lists.newArrayList();
        fromVCFFile(tumor, reference, rna, vcfFile, true, result::add);
        return result;
    }

    public void fromVCFFile(
            final String tumor, @Nullable final String reference, @Nullable final String rna,
            final String vcfFile, boolean useCheckReference, final Consumer<SmallVariant> consumer) throws IOException
    {
        VcfFileReader reader = new VcfFileReader(vcfFile);

        final VCFHeader header = reader.vcfHeader();

        if(!sampleInFile(tumor, header))
        {
            throw new IllegalArgumentException("Sample " + tumor + " not found in vcf file " + vcfFile);
        }

        if(useCheckReference && reference != null && !sampleInFile(reference, header))
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
            if(!mFilter.test(variantContext))
            {
                ++mFilteredCount;
                continue;
            }

            Optional<SmallVariant> varOptional = createVariant(tumor, reference, rna, variantContext);

            if(!varOptional.isPresent())
            {
                ++mFilteredCount;
                continue;
            }

            SmallVariant variant = varOptional.get();

            if(mDropDuplicates)
            {
                if(mLastNonFilteredVariant != null && variantsMatch(mLastNonFilteredVariant, variant))
                {
                    ++mFilteredCount;
                    continue;
                }

                mLastNonFilteredVariant = variant;
            }

            consumer.accept(variant);
            ++mCreatedCount;
        }
    }

    private static boolean variantsMatch(final SmallVariant first, final SmallVariant second)
    {
        return first.chromosome().equals(second.chromosome())
                && first.position() == second.position()
                && first.ref().equals(second.ref()) && first.alt().equals(second.alt());
    }

    public Optional<SmallVariant> createVariant(final String sample, final VariantContext context)
    {
        return createVariant(sample, null, null, context);
    }

    public Optional<SmallVariant> createVariant(
            final String sample, @Nullable final String reference, @Nullable final String rna, final VariantContext context)
    {
        final Genotype genotype = context.getGenotype(sample);

        if(mFilter.test(context) && AllelicDepth.containsAllelicDepth(genotype))
        {
            AllelicDepth tumorDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));

            Optional<AllelicDepth> referenceDepth = Optional.ofNullable(reference)
                    .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                    .filter(AllelicDepth::containsAllelicDepth)
                    .map(AllelicDepth::fromGenotype);

            if(tumorDepth.TotalReadCount > 0)
            {
                ImmutableSmallVariantImpl.Builder builder = createVariantBuilder(tumorDepth, context, reference, rna);

                return Optional.of(builder)
                        .map(x -> x.referenceDepth(referenceDepth.orElse(null)))
                        .map(ImmutableSmallVariantImpl.Builder::build);
            }
        }
        return Optional.empty();
    }

    public static ImmutableSmallVariantImpl.Builder createVariantBuilder(
            final AllelicDepth allelicDepth, final VariantContext context,
            @Nullable final String reference, @Nullable final String rna)
    {
        VariantContextDecorator decorator = new VariantContextDecorator(context);

        AllelicDepth rnaDepth = extractRnaDepth(context, rna);
        VariantImpact variantImpact = decorator.variantImpact();
        GenotypeStatus genotypeStatus = reference != null ? decorator.genotypeStatus(reference) : GenotypeStatus.UNKNOWN;

        PathogenicSummary pathogenicSummary = context.hasAttribute(CLNSIG) ? PathogenicSummaryFactory.fromContext(context) : null;

        ImmutableSmallVariantImpl.Builder builder = ImmutableSmallVariantImpl.builder()
                .chromosome(decorator.chromosome())
                .position(decorator.position())
                .type(decorator.type())
                .gene(variantImpact.GeneName)
                .ref(decorator.ref())
                .alt(decorator.alt())
                .allelicDepth(allelicDepth)
                .canonicalTranscript(variantImpact.CanonicalTranscript)
                .canonicalEffect(variantImpact.CanonicalEffect)
                .canonicalCodingEffect(variantImpact.CanonicalCodingEffect)
                .canonicalHgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .canonicalHgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .qual(decorator.qual())
                .mappability(decorator.mappability())
                .filter(decorator.filter())
                .genesAffected(variantImpact.GenesAffected)
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .otherReportedEffects(variantImpact.OtherReportableEffects)
                .worstCodingEffect(variantImpact.WorstCodingEffect)
                .hotspot(decorator.hotspot())
                .adjustedCopyNumber(decorator.adjustedCopyNumber())
                .adjustedVAF(decorator.adjustedVaf())
                .minorAlleleCopyNumber(decorator.minorAlleleCopyNumber())
                .variantCopyNumber(decorator.variantCopyNumber())
                .biallelic(decorator.biallelic())
                .reported(decorator.reported())
                .genotypeStatus(genotypeStatus)
                .germlineStatus(extractGermlineStatus(context))
                .trinucleotideContext(decorator.trinucleotideContext())
                .microhomology(decorator.microhomology())
                .repeatSequence(decorator.repeatSequence())
                .repeatCount(decorator.repeatCount())
                .tier(decorator.tier())
                .rnaDepth(rnaDepth)
                .reportableTranscripts(decorator.reportableTranscripts())
                .clinvarInfo(context.getAttributeAsString(CLNSIG, ""))
                .subclonalLikelihood(context.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0))
                .kataegis(context.getAttributeAsString(KATAEGIS_FLAG, Strings.EMPTY))
                .gnomadFrequency(context.getAttributeAsDouble(GNOMAD_FREQ, 0))
                .somaticLikelihood(SomaticLikelihood.valueOf(
                        context.getAttributeAsString(PANEL_SOMATIC_LIKELIHOOD, SomaticLikelihood.UNKNOWN.toString())))
                .pathogenic(pathogenicSummary != null ? pathogenicSummary.Status.isPathogenic() : false)
                .pathogenicity(pathogenicSummary != null ? pathogenicSummary.Status.toString() : "");

        if(context.hasAttribute(LOCAL_PHASE_SET))
        {
            builder.localPhaseSets(context.getAttributeAsIntList(LOCAL_PHASE_SET, 0));
        }

        return builder;
    }

    private static boolean sampleInFile(final String sample, final VCFHeader header)
    {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    private static GermlineStatus extractGermlineStatus(final VariantContext context)
    {
        try
        {
            return GermlineStatus.valueOf(context.getAttributeAsString(PURPLE_GERMLINE_INFO, GermlineStatus.UNKNOWN.name()));
        }
        catch(IllegalArgumentException ignored)
        {
            // protect against changes in types, which aren't expected to impact downstream processing
            return GermlineStatus.UNKNOWN;
        }
    }

    @Nullable
    private static AllelicDepth extractRnaDepth(VariantContext context, @Nullable String rna)
    {
        return Optional.ofNullable(context.getGenotype(rna))
                .filter(AllelicDepth::containsAllelicDepth)
                .map(AllelicDepth::fromGenotype)
                .orElse(null);
    }

    @Override
    public boolean test(final VariantContext variantContext)
    {
        return mFilter.test(variantContext);
    }
}
