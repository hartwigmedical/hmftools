package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.DEL_STR;

import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.lilac.evidence.AminoAcid;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;

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

        for(Map.Entry<Integer, Nucleotide> entry2 : frag2.nucleotidesByLoci().entrySet())
        {
            int locus2 = entry2.getKey();
            Nucleotide nuc2 = entry2.getValue();
            if(frag1.containsNucleotideLocus(locus2))
                continue;

            frag1.addNucleotide(nuc2);
            frag1.addReads(frag2);
        }

        return frag1;
    }

    public static Fragment copyNucleotideFragment(final Fragment fragment)
    {
        // ignores all state, just starts with original information
        Fragment newFragment = new Fragment(fragment.reads().get(0), fragment.readGene(), fragment.genes(), fragment.rawNucleotidesByLoci().values());

        for(int i = 1; i < fragment.reads().size(); ++i)
        {
            newFragment.addRead(fragment.reads().get(i));
        }

        return newFragment;
    }

    public static String formCodonAminoAcid(int locus, final SortedMap<Integer, Nucleotide> nucleotides)
    {
        SortedMap<Integer, Nucleotide> head = nucleotides.tailMap(locus * 3);
        List<String> bases = head.values().stream().limit(3).map(Nucleotide::bases).toList();
        String first = bases.get(0);
        String second = bases.get(1);
        String third = bases.get(2);

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
        if(fragment.nucleotidesByLoci().isEmpty() || fragment.rawNucleotidesByLoci().isEmpty())
        {
            LL_LOGGER.warn("{} {} has no bases", fragment.id(), fragment.readInfo());
            return false;
        }

        if(fragment.aminoAcidConversionCount() > 0)
        {
            int i = 0;
            for(AminoAcid aminoAcid : fragment.aminoAcidsByLoci().values())
            {
                if(aminoAcid.acid().isEmpty())
                {
                    LL_LOGGER.warn("{} {} empty amino-acid, index({})", fragment.id(), fragment.readInfo(), i);
                    return false;
                }

                i++;
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
}
