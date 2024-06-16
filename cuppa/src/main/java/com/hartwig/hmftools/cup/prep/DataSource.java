package com.hartwig.hmftools.cup.prep;

import com.hartwig.hmftools.cup.common.CategoryType;

public enum DataSource
{
    DNA,
    RNA;

    public static DataSource fromCategory(final CategoryType categoryType)
    {
        switch(categoryType)
        {
            case GENE_EXP:
            case ALT_SJ:
                return RNA;
            default:
                return DNA;
        }
    }
}
