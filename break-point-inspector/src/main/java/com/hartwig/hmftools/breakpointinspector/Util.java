package com.hartwig.hmftools.breakpointinspector;

import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.breakpointinspector.datamodel.HMFVariantType;
import com.hartwig.hmftools.breakpointinspector.datamodel.Range;

class HMFVariantContext {
    private static final int SHORT_VARIANT_LENGTH = 1000;

    String Id;
    Location MantaBP1;
    Range Uncertainty1;
    Location MantaBP2;
    Range Uncertainty2;
    HMFVariantType Type;
    boolean Imprecise;
    boolean BND;
    HashSet<String> Filter = Sets.newHashSet();
    String InsertSequence;
    String HomologySequence;

    int OrientationBP1 = 0;
    int OrientationBP2 = 0;

    HMFVariantContext(final String id, final Location bp1, final Location bp2, final HMFVariantType type, final boolean imprecise) {
        Id = id;
        MantaBP1 = bp1;
        MantaBP2 = bp2;
        Type = type;
        Imprecise = imprecise;
    }

    private boolean isShortDelete() {
        boolean b = Type == HMFVariantType.DEL;
        b &= MantaBP1.ReferenceIndex == MantaBP2.ReferenceIndex;
        b &= (MantaBP2.Position - MantaBP1.Position) < SHORT_VARIANT_LENGTH;
        return b;
    }

    private boolean isShortDuplicate() {
        boolean b = Type == HMFVariantType.DUP;
        b &= MantaBP1.ReferenceIndex == MantaBP2.ReferenceIndex;
        b &= (MantaBP2.Position - MantaBP1.Position) < SHORT_VARIANT_LENGTH;
        return b;
    }

    boolean isShortVariant() {
        return isShortDelete() || isShortDuplicate();
    }

    boolean isInsert() {
        return Type == HMFVariantType.INS;
    }
}

class Util {

    static List<String> prefixList(final List<String> list, final String prefix) {
        return list.stream().map(s -> prefix + s).collect(Collectors.toList());
    }

    static <U> List<String> toStrings(final List<U> list) {
        return list.stream().map(Objects::toString).collect(Collectors.toList());
    }

}