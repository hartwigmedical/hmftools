package com.hartwig.hmftools.pave;

public enum SpliceImpactType
{
    UNKNOWN,
    HOMOLOGY_SHIFT, // ref bases before or after an INDEL can maintain the ref splice acceptor or donor motif
    BASE_SHIFT, // an INDEL disrupts the splice acceptor or donor bases
    BASE_CHANGE, // a MNV changes the splice acceptor or donor bases
    REGION_DELETED, // a DEL deletes the entire splice region
    OUTSIDE_RANGE; // the variant has no impact on the splice region

    public boolean isDisruptive()
    {
        switch(this)
        {
            case BASE_CHANGE:
            case BASE_SHIFT:
            case REGION_DELETED:
            case HOMOLOGY_SHIFT:
                return true;

            case OUTSIDE_RANGE:
            default:
                return false;

        }
    }
}
