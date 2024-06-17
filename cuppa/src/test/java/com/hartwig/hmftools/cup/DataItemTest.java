package com.hartwig.hmftools.cup;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.DataItemMatrix;
import com.hartwig.hmftools.cup.prep.DataSource;
import com.hartwig.hmftools.cup.prep.ItemType;

public class DataItemTest
{
    @Test
    public void canSortIndex()
    {
        List<DataItem.Index> indexes = Arrays.asList(
                new DataItem.Index(DataSource.DNA, ItemType.SIGNATURE, "SIG_2_13"),
                new DataItem.Index(DataSource.DNA, ItemType.SIGNATURE, "SIG_1"),

                new DataItem.Index(DataSource.DNA, ItemType.GEN_POS, "X_0"),
                new DataItem.Index(DataSource.DNA, ItemType.GEN_POS, "10_0"),
                new DataItem.Index(DataSource.DNA, ItemType.GEN_POS, "2_0"),

                new DataItem.Index(DataSource.RNA, ItemType.ALT_SJ, "X;0;1"),
                new DataItem.Index(DataSource.RNA, ItemType.ALT_SJ, "10;0;10"),
                new DataItem.Index(DataSource.RNA, ItemType.ALT_SJ, "10;0;1"),
                new DataItem.Index(DataSource.RNA, ItemType.ALT_SJ, "2;0;1")
        );

        Collections.sort(indexes);

        String[] actualSortedKeys = indexes.stream().map(o -> o.Key).toArray(String[]::new);
        String[] expectedSortedKeys = {
                "SIG_1", "SIG_2_13",
                "2_0", "10_0", "X_0",
                "2;0;1", "10;0;1", "10;0;10", "X;0;1"
        };

        assertEquals(expectedSortedKeys, actualSortedKeys);
    }

    @Test
    public void canCreateDataItemMatrixWithUniqueKeys()
    {
        List<String> sampleIds = Arrays.asList("Sample1", "Sample2");
        ConcurrentHashMap<DataItem.Index, String[]> featureBySampleMatrix = new ConcurrentHashMap<>();

        List<List<DataItem>> dataItemsPerSample = Arrays.asList(
                List.of(
                        new DataItem(DataSource.DNA, ItemType.DRIVER, "BRAF.mut", "1"),
                        new DataItem(DataSource.DNA, ItemType.DRIVER, "TP53.mut", "0.5"),
                        new DataItem(DataSource.DNA, ItemType.DRIVER, "TERT.mut", "1")
                ),

                List.of(
                        new DataItem(DataSource.DNA, ItemType.DRIVER, "BRAF.mut", "1")
                )
        );

        int sampleIndex = 0;
        for(List<DataItem> dataItems : dataItemsPerSample)
        {
            for(DataItem dataItem : dataItems)
            {
                featureBySampleMatrix.computeIfAbsent(dataItem.Index, k -> new String[sampleIds.size()]);
                featureBySampleMatrix.get(dataItem.Index)[sampleIndex] = dataItem.Value;
            }
            sampleIndex++;
        }

        DataItemMatrix dataItemMatrix = new DataItemMatrix(sampleIds, featureBySampleMatrix);

        int nKeys = (int) dataItemMatrix.FeatureBySampleMatrix.keySet()
                .stream()
                .distinct()
                .count();

        assertEquals(3, nKeys);

        assertEquals(
                new String[] { "1", "1" },
                dataItemMatrix.get(new DataItem.Index(DataSource.DNA, ItemType.DRIVER, "BRAF.mut"))
        );

        assertEquals(
                new String[] { "0.5", null },
                dataItemMatrix.get(new DataItem.Index(DataSource.DNA, ItemType.DRIVER, "TP53.mut"))
        );
    }
}
