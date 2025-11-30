package com.hartwig.hmftools.orange.algo.purple;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS_DELIM;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantBuilderUtils;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.orange.algo.util.VariantTranscriptImpactCleaner;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFHeader;

public class PurpleVariantContextLoader
{
    private final CompoundFilter mFilter;

    public static PurpleVariantContextLoader withPassingOnlyFilter()
    {
        return new PurpleVariantContextLoader(new PassingVariantFilter());
    }

    public PurpleVariantContextLoader(final VariantContextFilter... filters)
    {
        mFilter = new CompoundFilter(true);
        mFilter.addAll(Arrays.asList(filters));
        mFilter.add(new HumanChromosomeFilter());
        mFilter.add(new NTFilter());
    }

    public List<PurpleVariantContext> fromVCFFile(
            final String tumor, @Nullable final String reference, @Nullable final String rna,
            final String vcfFile) throws IOException
    {
        List<PurpleVariantContext> result = new ArrayList<>();

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        final VCFHeader header = vcfFileReader.vcfHeader();

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

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(mFilter.test(variantContext))
            {
                PurpleVariantContext purpleVariantContext = createPurpleVariantContext(variantContext, tumor, reference, rna);
                result.add(purpleVariantContext);
            }
        }
        return result;
    }

    public PurpleVariantContext createPurpleVariantContext(
            final VariantContext variantContext, final String sample, @Nullable String reference, @Nullable String rna)
    {
        if(!AllelicDepth.containsAllelicDepth(variantContext.getGenotype(sample)))
        {
            throw new IllegalArgumentException(String.format(
                    "Variant could not be created because sample [%s] does not contain allelic depth",
                    sample));
        }

        final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(variantContext.getGenotype(sample));
        int readCount = tumorDepth.TotalReadCount;

        if(readCount < 0)
        {
            throw new IllegalArgumentException(String.format(
                    "Variant could not be created because tumor depth read count was negative! (value: %s)",
                    readCount));
        }

        return createPurpleVariantContext(variantContext, tumorDepth, reference, rna);
    }

    private PurpleVariantContext createPurpleVariantContext(
            final VariantContext variantContext, AllelicDepth tumorDepth, @Nullable String reference, @Nullable String rna)
    {
        VariantContextDecorator contextDecorator = new VariantContextDecorator(variantContext);

        final VariantImpact variantImpact = contextDecorator.variantImpact();

        final List<VariantTranscriptImpact> allImpacts = VariantTranscriptImpact.fromVariantContext(variantContext)
                .stream()
                .map(VariantTranscriptImpactCleaner::cleanFields)
                .collect(Collectors.toList());

        final List<VariantTranscriptImpact> otherImpacts = filterOutCanonicalImpact(allImpacts, variantImpact.CanonicalTranscript);

        return ImmutablePurpleVariantContext.builder()
                .variant(VariantBuilderUtils.createVariantBuilder(tumorDepth, variantContext, reference, rna).build())
                .otherImpacts(otherImpacts)
                .biallelicProbability(contextDecorator.biallelicProbability())
                .subclonalLikelihood(variantContext.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0))
                .build();
    }

    private static List<VariantTranscriptImpact> filterOutCanonicalImpact(
            final List<VariantTranscriptImpact> impacts, String canonicalTranscript)
    {
        return impacts.stream().filter(impact -> !impact.Transcript.equals(canonicalTranscript)).collect(Collectors.toList());
    }

    private static boolean sampleInFile(final String sample, final VCFHeader header)
    {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }
}
