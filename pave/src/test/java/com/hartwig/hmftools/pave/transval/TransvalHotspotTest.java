package com.hartwig.hmftools.pave.transval;

import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.Assert;
import org.junit.Test;

public class TransvalHotspotTest
{
    @Test
    public void compareToTest()
    {
        TransvalHotspot th0 = new TransvalHotspot("GGGC", "T", "chr5", 10_000_123);
        TransvalHotspot th1 = new TransvalHotspot("GGG", "ACC", "chr5", 10_000_123);
        TransvalHotspot th2 = new TransvalHotspot("GGG", "ACG", "chr5", 10_000_123);
        TransvalHotspot th3 = new TransvalHotspot("GGG", "ACG", "chr5", 10_000_124);
        TransvalHotspot th4 = new TransvalHotspot("GGGA", "ACCC", "chr5", 10_000_123);
        TransvalHotspot th5 = new TransvalHotspot("GGGA", "ACCC", "chr5", 10_000_124);
        TransvalHotspot th6 = new TransvalHotspot("GGGAA", "ACCC", "chr5", 10_000_124);
        TransvalHotspot th7 = new TransvalHotspot("GGGAAA", "ACCC", "chr5", 10_000_124);
        SortedSet<TransvalHotspot> sortedSet = new TreeSet<>();
        sortedSet.add(th7);
        sortedSet.add(th6);
        sortedSet.add(th5);
        sortedSet.add(th4);
        sortedSet.add(th3);
        sortedSet.add(th2);
        sortedSet.add(th1);
        sortedSet.add(th0);
        Iterator<TransvalHotspot> iterator = sortedSet.iterator();
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
