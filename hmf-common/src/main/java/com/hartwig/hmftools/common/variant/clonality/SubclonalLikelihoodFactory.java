package com.hartwig.hmftools.common.variant.clonality;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public class SubclonalLikelihoodFactory {

    private static final double MAX_PLOIDY = 2;

    public static List<SubclonalLikelihood> subclonalLikelihood(double binWidth, @NotNull final List<PeakModel> peakModel) {
        final List<SubclonalLikelihood> result = Lists.newArrayList();
        final Map<Integer, ModifiableSubclonalLikelihood> clonalityMap = Maps.newHashMap();

        for (final PeakModel entry : peakModel) {
            if (entry.isValid() && Doubles.lessThan(entry.peak(), MAX_PLOIDY)) {
                int key = (int) Math.round(entry.bucket() / binWidth);
                final ModifiableSubclonalLikelihood likelihood = clonalityMap.computeIfAbsent(key,
                        integer -> ModifiableSubclonalLikelihood.create()
                                .setBucket(entry.bucket())
                                .setClonalWeight(0)
                                .setSubclonalWeight(0));

                if (entry.isSubclonal()) {
                    likelihood.setSubclonalWeight(likelihood.subclonalWeight() + entry.bucketWeight());
                } else {
                    likelihood.setClonalWeight(likelihood.clonalWeight() + entry.bucketWeight());
                }
            }
        }

        for (final ModifiableSubclonalLikelihood entry : clonalityMap.values()) {
            double likelihood = entry.subclonalWeight() / (entry.subclonalWeight() + entry.clonalWeight());
            if (!Doubles.isZero(likelihood)) {
                result.add(entry.setLikelihood(likelihood));
            }
        }

        return result;
    }

}
