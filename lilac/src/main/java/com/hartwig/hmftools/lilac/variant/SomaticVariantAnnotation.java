package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.LilacConstants;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SomaticVariantAnnotation
{
    private final LilacConfig mConfig;

    // variant loci are grouped by the genes they nominally fall within
    private final Map<String,List<Integer>> mGeneVariantLoci;

    private final List<SomaticVariant> mSomaticVariants;

    private final Set<CodingEffect> UNKNOWN_CODING_EFFECT;
    private final Map<String,TranscriptData> mHlaTranscriptData;

    private final LociPosition mLociPositionFinder;

    public SomaticVariantAnnotation(
            final LilacConfig config, final Map<String,TranscriptData> transcriptData, final LociPosition lociPositionFinder)
    {
        mConfig = config;
        mHlaTranscriptData = transcriptData;
        UNKNOWN_CODING_EFFECT = Sets.newHashSet(CodingEffect.NONE, CodingEffect.UNDEFINED);

        mGeneVariantLoci = Maps.newHashMap();

        mLociPositionFinder = lociPositionFinder;

        mSomaticVariants = Lists.newArrayList();

        List<SomaticVariant> variants = loadSomaticVariants();

        for(SomaticVariant variant : variants)
        {
            int variantNucleotideLoci = mLociPositionFinder.calcNucelotideLocus(variant.Position);

            if(variantNucleotideLoci < 0)
                continue;

            mSomaticVariants.add(variant);

            int variantAminoAcidLoci = variantNucleotideLoci / 3;

            List<Integer> geneLoci = mGeneVariantLoci.get(variant.Gene);

            if(geneLoci == null)
            {
                geneLoci = Lists.newArrayList();
                mGeneVariantLoci.put(variant.Gene, geneLoci);
            }

            geneLoci.add(variantAminoAcidLoci);
        }
    }

    public List<SomaticVariant> getSomaticVariants() { return mSomaticVariants; }

    public final List<AlleleCoverage> assignAlleleCoverage(
            final SomaticVariant variant, final BamReader reader, final List<HlaSequenceLoci> winners)
    {
        List<Fragment> fragments = reader.readFromBam(variant);

        List<AlleleCoverage> coverages = Lists.newArrayList();

        if(fragments.isEmpty())
            return coverages;

        fragments.forEach(x -> x.qualityFilter(mConfig.MinBaseQual));
        fragments.forEach(x -> x.buildAminoAcids());

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

                    int index = fragment.getAminoAcidLoci().indexOf(locus);
                    String fragmentAA = "";

                    if(index >= 0)
                    {
                        fragmentAA = fragment.getAminoAcids().get(index);
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

        int maxWildcard = coverages.stream().mapToInt(x -> (int)x.WildCoverage).max().orElse(0);

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
            VCFFileReader fileReader = new VCFFileReader(new File(mConfig.SomaticVariantsFile), false);

            final CloseableIterator<VariantContext> variantIter = fileReader.iterator();

            while(variantIter.hasNext())
            {
                VariantContext variantContext = variantIter.next();

                if(variantContext.isFiltered() && !variantContext.getFilters().contains(SomaticVariantFactory.PASS_FILTER))
                    continue;

                VariantContextDecorator variant = new VariantContextDecorator(variantContext);

                if(inHlaCodingRegion(variant))
                {
                    variants.add(new SomaticVariant(
                            variant.gene(), variant.chromosome(), variant.position(), variant.ref(), variant.alt(),
                            variant.filter(), variant.canonicalCodingEffect(), variant.context()));
                }
            }

            LL_LOGGER.info("loaded {} HLA variants from file: {}", variants.size(), mConfig.SomaticVariantsFile);
            fileReader.close();
        }
        else if(!mConfig.CohortSomaticVariantsFile.isEmpty())
        {
            try
            {
                // load a cohort file - for now only retain the required sample's data
                final List<String> fileData = Files.readAllLines(new File(mConfig.SomaticVariantsFile).toPath());
                String header = fileData.get(0);
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);
                fileData.remove(0); // remove header
                int sampleIndex = fieldsIndexMap.get("SampleId");
                int geneIndex = fieldsIndexMap.get("Gene");
                int chrIndex = fieldsIndexMap.get("Chromosome");
                int posIndex = fieldsIndexMap.get("Position");
                int refIndex = fieldsIndexMap.get("Ref");
                int altIndex = fieldsIndexMap.get("Alt");
                int filterIndex = fieldsIndexMap.get("Filter");
                int effectIndex = fieldsIndexMap.get("CanonicalCodingEffect");

                for(final String line : fileData)
                {
                    String[] items = line.split(DELIM, -1);

                    String sampleId = items[sampleIndex];

                    if(!sampleId.equals(mConfig.Sample))
                        continue;

                    String gene = items[geneIndex];
                    CodingEffect codingEffect = CodingEffect.valueOf(items[effectIndex]);
                    String filter = items[filterIndex];

                    if(!HLA_GENES.contains(gene))
                        continue;

                    if(UNKNOWN_CODING_EFFECT.contains(codingEffect))
                        continue;

                    if(!filter.equals(SomaticVariantFactory.PASS_FILTER))
                        continue;

                    variants.add(new SomaticVariant(
                            gene, items[chrIndex], Integer.parseInt(items[posIndex]), items[refIndex], items[altIndex],
                            filter, codingEffect, null));
                }

                LL_LOGGER.info("loaded {} HLA variants from cohort file: {}", variants.size(), mConfig.CohortSomaticVariantsFile);
            }
            catch(IOException e)
            {
                LL_LOGGER.error("failed to read somatic variant file({}): {}", mConfig.SomaticVariantsFile, e.toString());
            }
        }

        return variants;
    }

}
