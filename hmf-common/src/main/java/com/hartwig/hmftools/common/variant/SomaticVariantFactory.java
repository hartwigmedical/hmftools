package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS_DELIM;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.StringJoiner;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantFactory implements VariantContextFilter
{
    @NotNull
    public static SomaticVariantFactory passOnlyInstance()
    {
        return new SomaticVariantFactory(new PassingVariantFilter());
    }

    public static final String MAPPABILITY_TAG = "MAPPABILITY";
    private static final String RECOVERED_FLAG = "RECOVERED";

    public static final String PASS_FILTER = "PASS";

    @NotNull
    private final CompoundFilter mFilter;
    private int mCreatedCount;
    private int mFilteredCount;

    public SomaticVariantFactory(final VariantContextFilter... filters)
    {
        mFilter = new CompoundFilter(true);
        mFilter.addAll(Arrays.asList(filters));
        mFilter.add(new HumanChromosomeFilter());
        mFilter.add(new NTFilter());
        mCreatedCount = 0;
        mFilteredCount = 0;
    }

    public int getCreatedCount()
    {
        return mCreatedCount;
    }

    public int getFilteredCount()
    {
        return mFilteredCount;
    }

    public List<SomaticVariant> fromVCFFile(final String tumor, final String vcfFile) throws IOException
    {
        return fromVCFFile(tumor, null, null, vcfFile);
    }

    @NotNull
    public List<SomaticVariant> fromVCFFile(final String tumor, @Nullable final String reference,
            final String vcfFile) throws IOException
    {
        final List<SomaticVariant> result = Lists.newArrayList();
        fromVCFFile(tumor, reference, null, vcfFile, true, result::add);
        return result;
    }

    public List<SomaticVariant> fromVCFFile(
            final String tumor, @Nullable final String reference, @Nullable final String rna, final String vcfFile) throws IOException
    {
        final List<SomaticVariant> result = Lists.newArrayList();
        fromVCFFile(tumor, reference, rna, vcfFile, true, result::add);
        return result;
    }

    public void fromVCFFile(
            final String tumor, @Nullable final String reference, @Nullable final String rna,
            final String vcfFile, boolean useCheckReference, Consumer<SomaticVariant> consumer) throws IOException
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

        for(VariantContext variant : reader.iterator())
        {
            if(mFilter.test(variant))
            {
                Optional<SomaticVariant> varOptional = createVariant(tumor, reference, rna, variant);

                if(varOptional.isPresent())
                {
                    varOptional.ifPresent(consumer);
                    ++mCreatedCount;
                }
                else
                {
                    ++mFilteredCount;
                }
            }
            else
            {
                ++mFilteredCount;
            }
        }
    }

    public Optional<SomaticVariant> createVariant(final String sample, final VariantContext context)
    {
        return createVariant(sample, null, null, context);
    }

    public Optional<SomaticVariant> createVariant(
            final String sample, @Nullable final String reference, @Nullable final String rna, final VariantContext context)
    {
        final Genotype genotype = context.getGenotype(sample);

        final VariantContextDecorator decorator = new VariantContextDecorator(context);
        final GenotypeStatus genotypeStatus = reference != null ? decorator.genotypeStatus(reference) : null;

        if(mFilter.test(context) && AllelicDepth.containsAllelicDepth(genotype))
        {
            AllelicDepth tumorDepth = AllelicDepth.fromGenotype(context.getGenotype(sample));

            Optional<AllelicDepth> referenceDepth = Optional.ofNullable(reference)
                    .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                    .filter(AllelicDepth::containsAllelicDepth)
                    .map(AllelicDepth::fromGenotype);

            Optional<AllelicDepth> rnaDepth = Optional.ofNullable(rna)
                    .flatMap(x -> Optional.ofNullable(context.getGenotype(x)))
                    .filter(AllelicDepth::containsAllelicDepth)
                    .map(AllelicDepth::fromGenotype);

            if(tumorDepth.TotalReadCount > 0)
            {
                ImmutableSomaticVariantImpl.Builder builder = createVariantBuilder(tumorDepth, context);
                builder.genotypeStatus(genotypeStatus != null ? genotypeStatus : GenotypeStatus.UNKNOWN);

                return Optional.of(builder)
                        .map(x -> x.rnaDepth(rnaDepth.orElse(null)))
                        .map(x -> x.referenceDepth(referenceDepth.orElse(null)))
                        .map(ImmutableSomaticVariantImpl.Builder::build);
            }
        }
        return Optional.empty();
    }

    private static ImmutableSomaticVariantImpl.Builder createVariantBuilder(final AllelicDepth allelicDepth, final VariantContext context)
    {
        final VariantContextDecorator decorator = new VariantContextDecorator(context);
        final VariantImpact variantImpact = decorator.variantImpact();

        List<String> reportableTranscripts = null;

        if(context.hasAttribute(REPORTABLE_TRANSCRIPTS))
        {
            String reportableTransStr = context.getAttributeAsString(REPORTABLE_TRANSCRIPTS, "");

            reportableTranscripts = Arrays.stream(reportableTransStr.split("\\" + REPORTABLE_TRANSCRIPTS_DELIM, -1))
                    .collect(Collectors.toList());
        }

        // to protect against changes in types, which aren't expected to impact downstream processing
        GermlineStatus germlineStatus = GermlineStatus.UNKNOWN;

        try
        {
            germlineStatus = GermlineStatus.valueOf(context.getAttributeAsString(PURPLE_GERMLINE_INFO, germlineStatus.toString()));
        }
        catch(Exception e) {}

        ImmutableSomaticVariantImpl.Builder builder = ImmutableSomaticVariantImpl.builder()
                .qual(decorator.qual())
                .type(decorator.type())
                .filter(decorator.filter())
                .chromosome(decorator.chromosome())
                .position(decorator.position())
                .ref(decorator.ref())
                .alt(decorator.alt())
                .allelicDepth(allelicDepth)
                .hotspot(decorator.hotspot())
                .minorAlleleCopyNumber(decorator.minorAlleleCopyNumber())
                .adjustedCopyNumber(decorator.adjustedCopyNumber())
                .adjustedVAF(decorator.adjustedVaf())
                .variantCopyNumber(decorator.variantCopyNumber())
                .mappability(decorator.mappability())
                .tier(decorator.tier())
                .trinucleotideContext(decorator.trinucleotideContext())
                .microhomology(decorator.microhomology())
                .repeatCount(decorator.repeatCount())
                .repeatSequence(decorator.repeatSequence())
                .reported(decorator.reported())
                .biallelic(decorator.biallelic())
                .gene(variantImpact.GeneName)
                .canonicalTranscript(variantImpact.CanonicalTranscript)
                .canonicalEffect(variantImpact.CanonicalEffect)
                .canonicalCodingEffect(variantImpact.CanonicalCodingEffect)
                .canonicalHgvsCodingImpact(variantImpact.CanonicalHgvsCoding)
                .canonicalHgvsProteinImpact(variantImpact.CanonicalHgvsProtein)
                .spliceRegion(variantImpact.CanonicalSpliceRegion)
                .otherReportedEffects(variantImpact.OtherReportableEffects)
                .reportableTranscripts(reportableTranscripts)
                .worstCodingEffect(variantImpact.WorstCodingEffect)
                .genesAffected(variantImpact.GenesAffected)
                .subclonalLikelihood(context.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0))
                .germlineStatus(germlineStatus)
                .kataegis(context.getAttributeAsString(KATAEGIS_FLAG, Strings.EMPTY))
                .recovered(context.getAttributeAsBoolean(RECOVERED_FLAG, false))
                .clinvarInfo(context.getAttributeAsString(CLNSIG, ""))
                .gnomadFrequency(context.getAttributeAsDouble(GNOMAD_FREQ, 0))
                .somaticLikelihood(SomaticLikelihood.valueOf(
                        context.getAttributeAsString(PANEL_SOMATIC_LIKELIHOOD, SomaticLikelihood.UNKNOWN.toString())));

        if(context.hasAttribute(LOCAL_PHASE_SET))
        {
            builder.localPhaseSets(context.getAttributeAsIntList(LOCAL_PHASE_SET, 0));
        }

        return builder;
    }

    public static final String LPS_DELIM = ";";

    public static List<Integer> localPhaseSetsStringToList(@Nullable final String localPhaseSetStr)
    {
        if(localPhaseSetStr == null || localPhaseSetStr.isEmpty())
        {
            return null;
        }

        List<Integer> localPhaseSets = Lists.newArrayList();
        Arrays.stream(localPhaseSetStr.split(LPS_DELIM)).forEach(x -> localPhaseSets.add(Integer.valueOf(x)));
        return localPhaseSets;
    }

    public static String localPhaseSetsStr(@Nullable final List<Integer> localPhaseSets)
    {
        if(localPhaseSets == null || localPhaseSets.isEmpty())
        {
            return "";
        }

        if(localPhaseSets.size() == 1)
        {
            return String.valueOf(localPhaseSets.get(0));
        }

        StringJoiner sj = new StringJoiner(LPS_DELIM);
        localPhaseSets.forEach(x -> sj.add(String.valueOf(x)));
        return sj.toString();
    }

    private static boolean sampleInFile(final String sample, final VCFHeader header)
    {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    @Override
    public boolean test(final VariantContext variantContext)
    {
        return mFilter.test(variantContext);
    }
}
