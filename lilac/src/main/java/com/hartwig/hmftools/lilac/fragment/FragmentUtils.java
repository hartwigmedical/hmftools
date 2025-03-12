package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;

public class FragmentUtils
{
    public static List<Fragment> mergeFragmentsById(final List<Fragment> fragments)
    {
        // merge paired reads but keep any single reads as well
        final List<Fragment> mergedFragments = Lists.newArrayList(); // by readId

        Map<String,List<Fragment>> readGroupFrags = Maps.newHashMap();

        for(Fragment fragment : fragments)
        {
            List<Fragment> idFrags = readGroupFrags.get(fragment.id());

            if(idFrags == null)
            {
                idFrags = Lists.newArrayList();
                readGroupFrags.put(fragment.id(), idFrags);
            }

            idFrags.add(fragment);
        }

        for(Map.Entry<String,List<Fragment>> readGroup : readGroupFrags.entrySet())
        {
            List<Fragment> idFrags = readGroup.getValue();

            if(idFrags.size() == 1)
            {
                mergedFragments.add(idFrags.get(0));
            }
            else if(idFrags.size() == 2)
            {
                mergedFragments.add(mergeFragments(idFrags.get(0), idFrags.get(1)));
            }
            else
            {
                Fragment mergedFrag = mergeFragments(idFrags.get(0), idFrags.get(1));

                for(int i = 2; i < idFrags.size(); ++i)
                {
                    mergedFrag = mergeFragments(mergedFrag, idFrags.get(i));
                }

                mergedFragments.add(mergedFrag);
            }
        }

        return mergedFragments;
    }

    public static List<Integer> expandIndices(final List<Integer> aminoAcidLoci)
    {
        List<Integer> nucleotideLoci = Lists.newArrayList();
        aminoAcidLoci.stream().forEach(x -> nucleotideLoci.addAll(Lists.newArrayList(3 * x, 3 * x + 1, 3 * x + 2)));
        return nucleotideLoci;
    }

    public static Fragment mergeFragments(final Fragment frag1, final Fragment frag2)
    {
        // merge frag 2 into 1 and ensure no repetition of loci
        frag2.genes().forEach(x -> frag1.genes().add(x));

        int index1 = 0;

        for(int index2 = 0; index2 < frag2.nucleotideLoci().size(); ++index2)
        {
            int locus2 = frag2.nucleotideLoci().get(index2);

            boolean add = true;

            for(; index1 < frag1.nucleotideLoci().size(); ++index1)
            {
                int locus1 = frag1.nucleotideLoci().get(index1);

                if(locus1 < locus2)
                    continue;

                if(locus1 == locus2)
                    add = false;

                break;
            }

            if(add)
            {
                try
                {
                    frag1.addNucleotideInfo(index1, locus2, frag2.nucleotides().get(index2), frag2.nucleotideQuality().get(index2));
                    frag1.addReads(frag2);
                    ++index1; // move past insertion point
                }
                catch(Exception e)
                {
                    LL_LOGGER.error("merging frag({}: {}) with frag({}: {})", frag1.id(), frag1.readInfo(), frag2.id(), frag2.readInfo());
                }
            }
        }

        return frag1;
    }

    public static Fragment copyNucleotideFragment(final Fragment fragment)
    {
        // ignores all state, just starts with original information
        Fragment newFragment = new Fragment(fragment.reads().get(0), fragment.readGene(), fragment.genes(),
                fragment.rawNucleotideLoci(), fragment.rawNucleotideQuality(), fragment.rawNucleotides());

        for(int i = 1; i < fragment.reads().size(); ++i)
        {
            newFragment.addRead(fragment.reads().get(i));
        }

        return newFragment;
    }

    public static String formCodonAminoAcid(int locus, final List<Integer> nucleotideLoci, final List<String> nucleotides)
    {
        int nucIndex = nucleotideLoci.indexOf(locus * 3);
        String first = nucleotides.get(nucIndex);
        String second = nucleotides.get(nucIndex + 1);
        String third = nucleotides.get(nucIndex + 2);

        if(first == DEL_STR && second == DEL_STR && third == DEL_STR)
            return DEL_STR;

        return Codons.aminoAcidFromBases(first + second + third);
    }

    public static List<Integer> calcAminoAcidIndices(int nucStartIndex, int nucEndIndex)
    {
        int start = nucStartIndex / 3 + (nucStartIndex % 3 == 0 ? 0 : 1);
        int end = (nucEndIndex + 1) / 3 - 1;

        if(end <= start)
        {
            List<Integer> range = Lists.newArrayList();
            range.add(start);
            range.add(end);
            return range;
        }

        return formRange(start, end);
    }

    public static boolean validateFragment(final Fragment fragment)
    {
        if(fragment.genes().isEmpty())
        {
            LL_LOGGER.warn("{} {} has no genes", fragment.id(), fragment.readInfo());
            return false;
        }

        // check loci are ordered and consistent with qualities and bases
        if(fragment.nucleotides().isEmpty() || fragment.rawNucleotides().isEmpty())
        {
            LL_LOGGER.warn("{} {} has no bases", fragment.id(), fragment.readInfo());
            return false;
        }

        if(fragment.nucleotides().size() != fragment.nucleotideQuality().size())
        {
            LL_LOGGER.warn("{} {} inconsistent bases loci({}) quals({})",
                    fragment.id(), fragment.readInfo(), fragment.nucleotides().size(), fragment.nucleotideQuality().size());
            return false;
        }

        if(!validateLociBases(fragment.id(), fragment.nucleotideLoci(), fragment.nucleotides()))
            return false;

        if(!validateLociBases(fragment.id(), fragment.rawNucleotideLoci(), fragment.rawNucleotides()))
            return false;

        if(fragment.aminoAcidConversionCount() > 0)
        {
            if(!validateLociBases(fragment.id(), fragment.aminoAcidLoci(), fragment.aminoAcids()))
                return false;

            for(int i = 0; i < fragment.aminoAcidLoci().size(); ++i)
            {
                if(fragment.aminoAcids().get(i).isEmpty())
                {
                    LL_LOGGER.warn("{} {} empty amino-acid, index({})", fragment.id(), fragment.readInfo(), i);
                    return false;
                }
            }

            if(fragment.aminoAcidConversionCount() > 1)
            {
                LL_LOGGER.warn("{} {} amino-acid conversion repeated({})",
                        fragment.id(), fragment.readInfo(), fragment.aminoAcidConversionCount());
                return false;
            }
        }

        return true;
    }

    public static boolean validateLociBases(final String id, final List<Integer> loci, final List<String> sequences)
    {
        if(sequences.size() != sequences.size())
        {
            LL_LOGGER.warn("frag({}) inconsistent loci({}) bases({})",
                    id, loci.size(), sequences.size());
            return false;
        }

        for(int i = 0; i < loci.size() - 1; ++i)
        {
            if(loci.get(i) >= loci.get(i + 1))
            {
                LL_LOGGER.warn("frag({}) loci out of order at index({})", id, i);
                return false;
            }
        }

        return true;
    }
}
