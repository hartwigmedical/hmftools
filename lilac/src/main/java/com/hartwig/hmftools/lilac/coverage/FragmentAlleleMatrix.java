package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class FragmentAlleleMatrix
{
    private final List<FragmentAlleles> mFragmentAlleles;
    private final List<HlaAllele> mAlleles;

    private final Map<HlaAllele,Integer> mAlleleIndexMap;

    private final int mAlleleCount;
    private final int mFragCount;

    private final boolean[][][] mMatrix;

    private static int FULL = 0;
    private static int PARTIAL = 1;
    private static int WILD = 2;
    private static int ITEMS = WILD + 1;

    public FragmentAlleleMatrix(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        mAlleles = alleles;
        mFragmentAlleles = fragmentAlleles;

        mAlleleIndexMap = Maps.newHashMap();

        mAlleleCount = alleles.size();
        mFragCount = fragmentAlleles.size();

        mMatrix = new boolean[mFragCount][mAlleleCount][ITEMS];

        buildAlleleFragmentMatrix();
    }

    private void buildAlleleFragmentMatrix()
    {
        for(int alleleIndex = 0; alleleIndex < mAlleles.size(); ++alleleIndex)
        {
            HlaAllele allele = mAlleles.get(alleleIndex);
            mAlleleIndexMap.put(allele, alleleIndex);
        }

        for(int fragIndex = 0; fragIndex < mFragmentAlleles.size(); ++fragIndex)
        {
            FragmentAlleles fragment = mFragmentAlleles.get(fragIndex);

            for(int type = FULL; type <= WILD; ++type)
            {
                List<HlaAllele> alleles = type == FULL ? fragment.getFull() : type == PARTIAL ? fragment.getPartial() : fragment.getWild();
                for(HlaAllele allele : alleles)
                {
                    Integer alleleIndex = mAlleleIndexMap.get(allele);

                    if(alleleIndex == null)
                        continue;

                    mMatrix[fragIndex][alleleIndex][type] = true;
                }
            }
        }
    }

    public List<HlaAlleleCoverage> create(final HlaComplex complex)
    {
        List<HlaAllele> alleles = complex.getAlleles();
        int alleleCount = alleles.size();

        int[] alleleIndices = new int[alleleCount];
        int[] uniqueCoverage = new int[alleleCount];
        double[] combinedCoverage = new double[alleleCount];
        double[] wildCoverage = new double[alleleCount];

        for(int i = 0; i < alleleCount; ++i)
        {
            Integer alleleIndex = mAlleleIndexMap.get(alleles.get(i));
            if(alleleIndex == null)
                return Lists.newArrayList();

            alleleIndices[i] = alleleIndex;
        }

        boolean[] full = new boolean[alleleCount];
        boolean[] partial = new boolean[alleleCount];
        boolean[] wild = new boolean[alleleCount];

        for(int fragIndex = 0; fragIndex < mFragCount; ++fragIndex)
        {
            int fullCount = 0;
            int fullAlleleIndex = -1;
            int partialCount = 0;
            int wildCount = 0;

            for(int i = 0; i < alleleCount; ++i)
            {
                full[i] = false;
                partial[i] = false;
                wild[i] = false;

                int alleleIndex = alleleIndices[i];

                if(mMatrix[fragIndex][alleleIndex][FULL])
                {
                    full[i] = true;
                    ++fullCount;
                    fullAlleleIndex = i;
                }

                if(mMatrix[fragIndex][alleleIndex][PARTIAL])
                {
                    partial[i] = true;
                    ++partialCount;
                }

                if(mMatrix[fragIndex][alleleIndex][WILD])
                {
                    wild[i] = true;
                    ++wildCount;
                }
            }

            if(fullCount == 1 && partialCount == 0)
            {
                ++uniqueCoverage[fullAlleleIndex];
            }
            else if(fullCount > 0 || partialCount > 0)
            {
                double contribution = 1.0 / (fullCount + partialCount + wildCount);

                for(int i = 0; i < alleleCount; ++i)
                {
                    if(full[i])
                        combinedCoverage[i] += contribution;

                    if(partial[i])
                        combinedCoverage[i] += contribution;

                    if(wild[i])
                        wildCoverage[i] += contribution;
                }
            }

            /*
            Set<HlaAllele> fullAlleles = fragment.getFull().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toSet());
            Set<HlaAllele> partialAlleles = fragment.getPartial().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toSet());
            Set<HlaAllele> wildAlleles = fragment.getWild().stream().map(x -> asAlleleGroup ? x.asAlleleGroup() : x).collect(Collectors.toSet());

            if(fullAlleles.size() == 1 && partialAlleles.isEmpty())
            {
                increment(uniqueCoverageMap, fullAlleles.iterator().next(), 1);
            }
            else
            {
                double contribution = 1.0 / (fullAlleles.size() + partialAlleles.size() + wildAlleles.size());
                fullAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                partialAlleles.forEach(x -> increment(combinedCoverageMap, x, contribution));
                wildAlleles.forEach(x -> increment(wildCoverageMap, x, contribution));
            }
             */
        }


        List<HlaAlleleCoverage> alleleCoverages = Lists.newArrayListWithExpectedSize(alleleCount);

        for(int i = 0; i < alleleCount; ++i)
        {
            alleleCoverages.add(new HlaAlleleCoverage(alleles.get(i), uniqueCoverage[i], combinedCoverage[i], wildCoverage[i]));
        }

        // no known reason to sort
        // Collections.sort(results, Collections.reverseOrder());

        return alleleCoverages;
    }
}
