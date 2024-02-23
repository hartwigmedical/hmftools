package com.hartwig.hmftools.cup.prep;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

public class DataItemTest
{
    @Test
    public void canSortIndex()
    {
        String DUMMY_VALUE = "";

        List<DataItem> dataItems = Arrays.asList(
            new DataItem(DataSource.DNA, ItemType.SIGNATURE, "SIG_2_13", DUMMY_VALUE),
            new DataItem(DataSource.DNA, ItemType.SIGNATURE, "SIG_1", DUMMY_VALUE),
            new DataItem(DataSource.DNA, ItemType.GEN_POS, "X_0", DUMMY_VALUE),
            new DataItem(DataSource.DNA, ItemType.GEN_POS, "10_0", DUMMY_VALUE),
            new DataItem(DataSource.DNA, ItemType.GEN_POS, "2_0", DUMMY_VALUE)
        );

        List<DataItem.Index> indexes = new ArrayList<>();
        for(DataItem dataItem : dataItems)
            indexes.add(dataItem.Index);

        Collections.sort(indexes, new DataItem.IndexComparator());

        String[] actualSortedKeys = indexes.stream().map(o -> o.Key).toArray(String[]::new);
        String[] expectedSortedKeys = {"SIG_1", "SIG_2_13", "2_0", "10_0", "X_0"};

        assertEquals(expectedSortedKeys, actualSortedKeys);
    }
}
