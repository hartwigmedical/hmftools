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
    private static int WILD = 1;
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
                List<HlaAllele> alleles = type == FULL ? fragment.getFull() : fragment.getWild();
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

    public List<AlleleCoverage> create(final HlaComplex complex)
    {
        List<HlaAllele> alleles = complex.Alleles;
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
        boolean[] wild = new boolean[alleleCount];

        for(int fragIndex = 0; fragIndex < mFragCount; ++fragIndex)
        {
            int fullCount = 0;
            int fullAlleleIndex = -1;
            int wildCount = 0;

            for(int i = 0; i < alleleCount; ++i)
            {
                full[i] = false;
                wild[i] = false;

                int alleleIndex = alleleIndices[i];

                if(mMatrix[fragIndex][alleleIndex][FULL])
                {
                    full[i] = true;
                    ++fullCount;
                    fullAlleleIndex = i;
                }

                if(mMatrix[fragIndex][alleleIndex][WILD])
                {
                    wild[i] = true;
                    ++wildCount;
                }
            }

            if(fullCount == 1 && wildCount == 0)
            {
                ++uniqueCoverage[fullAlleleIndex];
            }
            else if(fullCount > 0 || wildCount > 0)
            {
                double contribution = 1.0 / (fullCount + wildCount);

                for(int i = 0; i < alleleCount; ++i)
                {
                    if(full[i])
                        combinedCoverage[i] += contribution;

                    if(wild[i])
                        wildCoverage[i] += contribution;
                }
            }
        }

        List<AlleleCoverage> alleleCoverages = Lists.newArrayListWithExpectedSize(alleleCount);

        for(int i = 0; i < alleleCount; ++i)
        {
            alleleCoverages.add(new AlleleCoverage(alleles.get(i), uniqueCoverage[i], combinedCoverage[i], wildCoverage[i]));
        }

        return alleleCoverages;
    }
}
