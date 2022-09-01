package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

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

    private final SupportType[][] mMatrix;

    private static int FULL = 0;
    private static int WILD = 1;
    private static int ITEMS = WILD + 1;

    private enum SupportType
    {
        FULL,
        WILD,
        BOTH,
        NONE;
    }

    public FragmentAlleleMatrix(final List<FragmentAlleles> fragmentAlleles, final List<HlaAllele> alleles)
    {
        mAlleles = alleles;
        mFragmentAlleles = fragmentAlleles;

        mAlleleIndexMap = Maps.newHashMap();

        mAlleleCount = alleles.size();
        mFragCount = fragmentAlleles.size();

        mMatrix = new SupportType[mFragCount][mAlleleCount];

        for(int f = 0; f < mFragCount; ++f)
        {
            for(int a = 0; a < mAlleleCount; ++a)
            {
                mMatrix[f][a] = SupportType.NONE;
            }
        }

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

            for(HlaAllele allele : fragment.getFull())
            {
                Integer alleleIndex = mAlleleIndexMap.get(allele);

                if(alleleIndex == null)
                    continue;

                mMatrix[fragIndex][alleleIndex] = SupportType.FULL;
            }

            for(HlaAllele allele : fragment.getWild())
            {
                Integer alleleIndex = mAlleleIndexMap.get(allele);

                if(alleleIndex == null)
                    continue;

                if(mMatrix[fragIndex][alleleIndex] == SupportType.FULL)
                    mMatrix[fragIndex][alleleIndex] = SupportType.BOTH;
                else
                    mMatrix[fragIndex][alleleIndex] = SupportType.WILD;
            }
        }
    }

    private class AlleleCounts
    {
        public final int AlleleIndex;

        public boolean Full = false;
        public boolean Wild = false;
        public int UniqueCoverage = 0;
        public double CombinedCoverage = 0;
        public double WildCoverage = 0;

        public AlleleCounts(int alleleIndex)
        {
            AlleleIndex = alleleIndex;
        }
    }

    public List<AlleleCoverage> create(final HlaComplex complex)
    {
        List<HlaAllele> alleles = complex.Alleles;
        int alleleCount = alleles.size();

        AlleleCounts[] alleleCounts = new AlleleCounts[alleleCount];

        for(int i = 0; i < alleleCount; ++i)
        {
            Integer alleleIndex = mAlleleIndexMap.get(alleles.get(i));
            if(alleleIndex == null)
                return Lists.newArrayList();

            alleleCounts[i] = new AlleleCounts(alleleIndex);
        }

        for(int fragIndex = 0; fragIndex < mFragCount; ++fragIndex)
        {
            int fullCount = 0;
            int fullAlleleIndex = -1;
            int wildCount = 0;

            for(int i = 0; i < alleleCount; ++i)
            {
                alleleCounts[i].Full = false;
                alleleCounts[i].Wild = false;

                int alleleIndex = alleleCounts[i].AlleleIndex;

                SupportType supportType = mMatrix[fragIndex][alleleIndex];

                if(supportType == SupportType.FULL || supportType == SupportType.BOTH)
                {
                    alleleCounts[i].Full = true;
                    ++fullCount;
                    fullAlleleIndex = i;
                }

                if(supportType == SupportType.WILD || supportType == SupportType.BOTH)
                {
                    alleleCounts[i].Wild = true;
                    ++wildCount;
                }
            }

            if(fullCount == 1 && wildCount == 0)
            {
                ++alleleCounts[fullAlleleIndex].UniqueCoverage;
            }
            else if(fullCount > 0 || wildCount > 0)
            {
                double contribution = 1.0 / (fullCount + wildCount);

                for(int i = 0; i < alleleCount; ++i)
                {
                    if(alleleCounts[i].Full)
                    {
                        alleleCounts[i].CombinedCoverage += contribution;
                    }

                    if(alleleCounts[i].Wild)
                    {
                        alleleCounts[i].WildCoverage += contribution;
                    }
                }
            }
        }

        List<AlleleCoverage> alleleCoverages = Lists.newArrayListWithExpectedSize(alleleCount);

        for(int i = 0; i < alleleCount; ++i)
        {
            alleleCoverages.add(new AlleleCoverage(
                    alleles.get(i), alleleCounts[i].UniqueCoverage, alleleCounts[i].CombinedCoverage, alleleCounts[i].WildCoverage));
        }

        return alleleCoverages;
    }
}
