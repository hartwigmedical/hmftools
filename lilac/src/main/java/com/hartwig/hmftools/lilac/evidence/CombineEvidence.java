package com.hartwig.hmftools.lilac.evidence;

import static java.lang.Math.min;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class CombineEvidence
{
    public static PhasedEvidence combine(final PhasedEvidence left, final PhasedEvidence common, final PhasedEvidence right)
    {
        PhasedEvidence leftWithCommon = combineOverlapping(left, common);
        PhasedEvidence result = combineOverlapping(leftWithCommon, right);
        return result;
    }

    public static boolean canCombine(final PhasedEvidence left, final PhasedEvidence common, final PhasedEvidence right)
    {
        return canCombine(left, common) && canCombine(common, right);
    }

    public static boolean canCombine(final PhasedEvidence left, final PhasedEvidence right)
    {
        List<Integer> indexIntersection = left.getAminoAcidLoci().stream()
                .filter(x -> right.getAminoAcidLoci().contains(x)).collect(Collectors.toList());

        Collections.sort(indexIntersection);

        if (indexIntersection.isEmpty())
            return false;

        List<Integer> uniqueToLeft = left.getAminoAcidLoci().stream().filter(x -> !indexIntersection.contains(x)).collect(Collectors.toList());
        List<Integer> uniqueToRight = right.getAminoAcidLoci().stream().filter(x -> !indexIntersection.contains(x)).collect(Collectors.toList());

        int rightMin = uniqueToRight.stream().mapToInt(x -> x).min().orElse(0);
        int leftMin = uniqueToLeft.stream().mapToInt(x -> x).min().orElse(0);
        if (!uniqueToRight.isEmpty() && uniqueToLeft.stream().anyMatch(x -> x > rightMin))
            return false;

        if (!uniqueToLeft.isEmpty() && uniqueToRight.stream().anyMatch(x -> x < leftMin))
            return false;

        Set<String> leftCommon = left.getEvidence().keySet().stream()
                .map(x -> x.substring(left.getAminoAcidLoci().size() - indexIntersection.size())).collect(Collectors.toSet());

        Set<String> rightCommon = right.getEvidence().keySet().stream()
                .map(x -> x.substring(0, indexIntersection.size())).collect(Collectors.toSet());

        return leftCommon.size() == rightCommon.size() && leftCommon.stream().allMatch(x -> rightCommon.contains(x));
    }

    public static PhasedEvidence combineOverlapping(final PhasedEvidence left, final PhasedEvidence right)
    {
        List<Integer> indexUnion = left.getAminoAcidLoci().stream().collect(Collectors.toList());
        right.getAminoAcidLoci().stream().filter(x -> !left.getAminoAcidLoci().contains(x)).forEach(x -> indexUnion.add(x));
        Collections.sort(indexUnion);

        List<Integer> indexIntersection = left.getAminoAcidLoci().stream()
                .filter(x -> right.getAminoAcidLoci().contains(x)).collect(Collectors.toList());

        Collections.sort(indexIntersection);

        List<Integer> uniqueToLeft = left.getAminoAcidLoci().stream().filter(x -> !indexIntersection.contains(x)).collect(Collectors.toList());

        Map<String,Integer> results = Maps.newHashMap();

        for(Map.Entry<String,Integer> leftEntry : left.getEvidence().entrySet())
        {

            String leftUnique = leftEntry.getKey().substring(0, uniqueToLeft.size());
            String leftCommon = leftEntry.getKey().substring(uniqueToLeft.size());

            for(Map.Entry<String, Integer> rightEntry : right.getEvidence().entrySet())
            {
                String rightCommon = rightEntry.getKey().substring(0, indexIntersection.size());

                if(leftCommon.equals(rightCommon))
                {
                    String combined = leftUnique + rightEntry.getKey();
                    results.put(combined, min(leftEntry.getValue(), rightEntry.getValue()));
                }
            }
        }

        return new PhasedEvidence(indexUnion, results);
    }
}
