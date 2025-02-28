package com.hartwig.hmftools.pave.transval;

import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.Assert;
import org.junit.Test;

public class BaseSequenceChangeTest
{
    @Test
    public void compareToTest()
    {
        BaseSequenceChange th0 = new BaseSequenceChange("GGGC", "T", "chr5", 10_000_123);
        BaseSequenceChange th1 = new BaseSequenceChange("GGG", "ACC", "chr5", 10_000_123);
        BaseSequenceChange th2 = new BaseSequenceChange("GGG", "ACG", "chr5", 10_000_123);
        BaseSequenceChange th3 = new BaseSequenceChange("GGG", "ACG", "chr5", 10_000_124);
        BaseSequenceChange th4 = new BaseSequenceChange("GGGA", "ACCC", "chr5", 10_000_123);
        BaseSequenceChange th5 = new BaseSequenceChange("GGGA", "ACCC", "chr5", 10_000_124);
        BaseSequenceChange th6 = new BaseSequenceChange("GGGAA", "ACCC", "chr5", 10_000_124);
        BaseSequenceChange th7 = new BaseSequenceChange("GGGAAA", "ACCC", "chr5", 10_000_124);
        SortedSet<BaseSequenceChange> sortedSet = new TreeSet<>();
        sortedSet.add(th7);
        sortedSet.add(th6);
        sortedSet.add(th5);
        sortedSet.add(th4);
        sortedSet.add(th3);
        sortedSet.add(th2);
        sortedSet.add(th1);
        sortedSet.add(th0);
        Iterator<BaseSequenceChange> iterator = sortedSet.iterator();
        Assert.assertEquals(th0, iterator.next());
        Assert.assertEquals(th1, iterator.next());
        Assert.assertEquals(th2, iterator.next());
        Assert.assertEquals(th3, iterator.next());
        Assert.assertEquals(th4, iterator.next());
        Assert.assertEquals(th5, iterator.next());
        Assert.assertEquals(th6, iterator.next());
        Assert.assertEquals(th7, iterator.next());
    }
}
