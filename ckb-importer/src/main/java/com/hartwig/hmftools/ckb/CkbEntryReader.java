package com.hartwig.hmftools.ckb;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.CkbJsonReader;

import org.jetbrains.annotations.NotNull;

public final class CkbEntryReader {

    private CkbEntryReader() {
    }

    @NotNull
    public static List<CkbEntry> read(@NotNull String ckbDir) throws IOException {
        CkbJsonDatabase ckbJsonDatabase = CkbJsonReader.read(ckbDir);
        return JsonDatabaseToCkbEntryConverter.convert(ckbJsonDatabase);
    }
}
