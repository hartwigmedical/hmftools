package com.hartwig.hmftools.orange.report.tables;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SignatureAllocationTableTest
{
    @Test
    public void canSortSignatureAllocations()
    {
        SignatureAllocation allocation1 = create("Sig1", 120D);
        SignatureAllocation allocation2 = create("Sig2", 200D);
        SignatureAllocation allocation3 = create("Sig10", 100D);
        SignatureAllocation allocation4 = create(SignatureAllocationTable.MISALLOC_SIGNATURE, 50D);

        List<SignatureAllocation> allocations = Lists.newArrayList(allocation2, allocation4, allocation3, allocation1);
        List<SignatureAllocation> sorted = SignatureAllocationTable.sort(allocations);

        assertEquals(4, sorted.size());
        assertEquals(allocation2, sorted.get(0));
        assertEquals(allocation1, sorted.get(1));
        assertEquals(allocation3, sorted.get(2));
        assertEquals(allocation4, sorted.get(3));
    }

    @NotNull
    private static SignatureAllocation create(@NotNull String signature, double allocation)
    {
        return ImmutableSignatureAllocation.builder().signature(signature).allocation(allocation).percent(0D).etiology("test").build();
    }
}