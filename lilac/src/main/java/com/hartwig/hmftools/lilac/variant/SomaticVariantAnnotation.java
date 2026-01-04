package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.calcNucelotideLocus;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.LilacConstants;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.evidence.AminoAcid;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantAnnotation
{
    private final LilacConfig mConfig;

    // variant loci are grouped by the genes they nominally fall within
    private final Map<HlaGene, List<Integer>> mGeneVariantLoci;

    private final List<SomaticVariant> mSomaticVariants;

    private final Map<HlaGene, TranscriptData> mHlaTranscriptData;

    public SomaticVariantAnnotation(final LilacConfig config, final Map<HlaGene, TranscriptData> transcriptData)
    {
        mConfig = config;
        mHlaTranscriptData = transcriptData;

        mGeneVariantLoci = Maps.newHashMap();

        List<TranscriptData> transcripts = transcriptData.values().stream().toList();

        mSomaticVariants = Lists.newArrayList();

        List<SomaticVariant> variants = loadSomaticVariants();

        for(SomaticVariant variant : variants)
        {
            int variantNucleotideLoci = calcNucelotideLocus(transcripts, variant.Position);

            if(variantNucleotideLoci < 0)
                continue;

            mSomaticVariants.add(variant);

            int variantAminoAcidLoci = variantNucleotideLoci / 3;

            mGeneVariantLoci.computeIfAbsent(variant.Gene, k -> Lists.newArrayList());
            List<Integer> geneLoci = mGeneVariantLoci.get(variant.Gene);
            geneLoci.add(variantAminoAcidLoci);
        }
    }

    public List<SomaticVariant> getSomaticVariants() { return mSomaticVariants; }

    public List<AlleleCoverage> assignAlleleCoverage(
            final SomaticVariant variant, final Collection<Fragment> fragments, final Iterable<HlaSequenceLoci> winners)
    {
        List<AlleleCoverage> coverages = Lists.newArrayList();

        if(fragments.isEmpty())
            return coverages;

        fragments.forEach(Fragment::removeLowQualBases);
        fragments.forEach(Fragment::buildAminoAcids);

        // take all variants in this gene together, and then ignore those loci
        List<Integer> variantLoci = mGeneVariantLoci.get(variant.Gene);

        for(HlaSequenceLoci sequenceLoci : winners)
        {
            if(variantLoci == null || variantLoci.isEmpty())
            {
                LL_LOGGER.error("no sequences found for variant({})", variant.toString());
                continue;
            }

            int supportCount = 0;

            for(Fragment fragment : fragments)
            {
                boolean matches = true;
                int matchCount = 0;
                for(int locus = fragment.minAminoAcidLocus(); locus <= fragment.maxAminoAcidLocus(); ++locus)
                {
                    if(locus >= sequenceLoci.length())
                        break;

                    if(variantLoci.contains(locus))
                        continue;

                    AminoAcid aminoAcid = fragment.aminoAcidsByLoci().get(locus);
                    String fragmentAA;

                    if(aminoAcid != null)
                    {
                        fragmentAA = aminoAcid.acid();
                    }
                    else
                    {
                        fragmentAA = fragment.getLowQualAminoAcid(locus);

                        if(fragmentAA.isEmpty())
                            continue;
                    }

                    String sequence = sequenceLoci.sequence(locus);

                    if(!fragmentAA.equals(sequence) && !sequence.equals(WILD_STR))
                    {
                        matches = false;
                        break;
                    }

                    ++matchCount;
                }

                if(matches && matchCount > 0)
                {
                    //LL_LOGGER.debug("allele({}) supported by fragment({} {}) matchedAAs({})",
                    //        seq.Allele, fragment.id(), fragment.readInfo(), matchCount);
                    ++supportCount;
                }
            }

            if(supportCount > 0)
            {
                if(sequenceLoci.hasWildcards())
                    coverages.add(new AlleleCoverage(sequenceLoci.Allele, 0, 0, supportCount));
                else
                    coverages.add(new AlleleCoverage(sequenceLoci.Allele, supportCount, 0, 0));
            }
        }

        // ignore wildcard alleles if any other have support from at least half the fragment
        int maxUnique = coverages.stream().mapToInt(x -> x.UniqueCoverage).max().orElse(0);

        if(maxUnique > 0 && maxUnique >= fragments.size() / 2)
            return coverages.stream().filter(x -> x.UniqueCoverage == maxUnique).collect(Collectors.toList());

        int maxWildcard = coverages.stream().mapToInt(x -> (int) x.WildCoverage).max().orElse(0);

        if(maxWildcard > 0)
            return coverages.stream().filter(x -> (int)x.WildCoverage == maxWildcard).collect(Collectors.toList());
        else
            return Lists.newArrayList();
    }

    private static class VariantCoverageSorter implements Comparator<AlleleCoverage>
    {
        // sorts by total coverage descending
        public int compare(final AlleleCoverage first, final AlleleCoverage second)
        {
            if(first.UniqueCoverage != second.UniqueCoverage)
                return first.UniqueCoverage < second.UniqueCoverage ? 1 : -1;

            if(first.WildCoverage != second.WildCoverage)
                return first.WildCoverage < second.WildCoverage ? 1 : -1;

            return 0;
        }
    }

    private boolean inHlaCodingRegion(final VariantContextDecorator variant)
    {
        if(!variant.chromosome().equals(LilacConstants.HLA_CHR))
            return false;

        int posStart = variant.position();
        int posEnd = posStart + variant.ref().length() - 1;

        for(TranscriptData transData : mHlaTranscriptData.values())
        {
            if(!positionsOverlap(posStart, posEnd, transData.CodingStart, transData.CodingEnd))
                continue;

            // check that the variant covers any exon
            if(transData.exons().stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.Start, x.End)))
                return true;

            // otherwise any splice variant
            if(variant.variantImpact() != null && variant.variantImpact().CanonicalCodingEffect == CodingEffect.SPLICE)
                return true;
        }

        return false;
    }

    private List<SomaticVariant> loadSomaticVariants()
    {
        List<SomaticVariant> variants = Lists.newArrayList();

        if(mConfig == null)
            return variants;

        if(!mConfig.SomaticVariantsFile.isEmpty())
        {
            VcfFileReader vcfReader = new VcfFileReader(mConfig.SomaticVariantsFile);

            for(VariantContext variantContext : vcfReader.iterator())
            {
                if(variantContext.isFiltered() && !variantContext.getFilters().contains(PASS_FILTER))
                    continue;

                VariantContextDecorator variant = new VariantContextDecorator(variantContext);

                if(inHlaCodingRegion(variant))
                {
                    variants.add(new SomaticVariant(
                            HlaGene.fromString(variant.gene()), variant.chromosome(), variant.position(), variant.ref(), variant.alt(),
                            variant.filter(), variant.canonicalCodingEffect(), variant.context()));
                }
            }

            LL_LOGGER.info("loaded {} HLA variants from file: {}", variants.size(), mConfig.SomaticVariantsFile);
            vcfReader.close();
        }

        return variants;
    }
}
