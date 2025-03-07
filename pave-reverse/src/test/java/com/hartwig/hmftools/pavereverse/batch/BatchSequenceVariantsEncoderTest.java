package com.hartwig.hmftools.pavereverse.batch;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.junit.Test;

public class BatchSequenceVariantsEncoderTest
{
    @Test
    public void columns()
    {
        List<String> expected = Lists.newArrayList("GENE", "HGVS_PROTEIN", "CHROMOSOME", "TRANSCRIPT", "VARIANTS");
        assertEquals(VariantsEncoder.columns(), expected);
    }

    @Test
    public void encode()
    {
        BaseSequenceChange change1 = new BaseSequenceChange("A", "AA", "chr2", 12345);
        BaseSequenceChange change2 = new BaseSequenceChange("T", "TT", "chr2", 12347);
        Set<BaseSequenceChange> changes = Sets.newHashSet();
        changes.add(change1);
        changes.add(change2);
        String str = VariantsEncoder.encodeChanges(changes);
        assertEquals(changes, parseChanges(str, "chr2"));
    }

    static Set<BaseSequenceChange> parseChanges(String str, String chr)
    {
        Set<BaseSequenceChange> changes = Sets.newHashSet();
        String[] split = str.split(";");
        Arrays.stream(split).forEach(s -> changes.add(parseChange(s, chr)));
        return changes;
    }

    static BaseSequenceChange parseChange(String s, String chr)
    {
        String[] parts = s.split(",");
        return new BaseSequenceChange(parts[0], parts[1], chr, Integer.parseInt(parts[2]));
    }
}
