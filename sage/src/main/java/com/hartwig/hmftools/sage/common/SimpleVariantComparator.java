package com.hartwig.hmftools.sage.common;

import java.util.Comparator;

public class SimpleVariantComparator implements Comparator<SimpleVariant>
{
    @Override
    public int compare(final SimpleVariant o1, final SimpleVariant o2)
    {
        int standardCompare = o1.compareTo(o2);

        if(standardCompare != 0)
            return standardCompare;

        int o1Length = Math.max(o1.ref().length(), o1.alt().length());
        int o2Length = Math.max(o2.ref().length(), o2.alt().length());
        int lengthCompare = Integer.compare(o1Length, o2Length);

        if(lengthCompare != 0)
            return lengthCompare;

        int refCompare = o1.ref().compareTo(o2.ref());

        if(refCompare != 0)
            return refCompare;

        return o1.alt().compareTo(o2.alt());
    }
}
