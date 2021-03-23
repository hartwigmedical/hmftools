package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;

import org.jetbrains.annotations.NotNull;

public class CKBFilter {

    @NotNull
    public List<CkbEntry> run(@NotNull List<CkbEntry> entries) {
        List<CkbEntry> filteredCKBEntries = Lists.newArrayList();
        return filteredCKBEntries;
    }

    public void reportUnusedFilterEntries() {

    }
}
