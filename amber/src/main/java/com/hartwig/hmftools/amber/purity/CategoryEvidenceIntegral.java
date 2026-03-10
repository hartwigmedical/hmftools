package com.hartwig.hmftools.amber.purity;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.Set;

public class CategoryEvidenceIntegral<T extends Comparable<T>>
{
    final List<CategoryEvidence<T>> mCategoryEvidence;

    public CategoryEvidenceIntegral(Set<CategoryEvidence<T>> categoryEvidence)
    {
        mCategoryEvidence = categoryEvidence.stream().sorted().collect(toList());
    }

    public int totalHits()
    {
        return mCategoryEvidence.stream().mapToInt(CategoryEvidence::evidencePoints).sum();
    }

    public int totalPoints()
    {
        return mCategoryEvidence.stream().map(CategoryEvidence::totalPoints).reduce(0, Integer::sum);
    }

    public double value()
    {
        double heightSoFar = 0.0;
        double area = 0.0;
        for(CategoryEvidence<T> categoryEvidence : mCategoryEvidence)
        {
            double rectangleSize = categoryEvidence.totalPoints() * heightSoFar;
            area += rectangleSize;
            double triangleSize = categoryEvidence.evidencePoints() * categoryEvidence.totalPoints() / 2.0;
            area += triangleSize;
            heightSoFar += categoryEvidence.evidencePoints();
        }
        return area;
    }
}
