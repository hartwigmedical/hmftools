package com.hartwig.hmftools.ecrfanalyser.reader;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ItemDef {
    @NotNull
    private final String OID;
    @NotNull
    private final String name;
    @Nullable
    private final String codeListOID;

    ItemDef(@NotNull final String OID, @NotNull final String name, @Nullable final String codeListOID) {
        this.OID = OID;
        this.name = name;
        this.codeListOID = codeListOID;
    }

    @NotNull
    String OID() {
        return OID;
    }

    @NotNull
    String name() {
        return name;
    }

    @Nullable
    String codeListOID() {
        return codeListOID;
    }

    @Override
    public String toString() {
        return "ItemDef{" + "OID='" + OID + '\'' + ", name='" + name + '\'' + ", codeListOID='" + codeListOID + '\''
                + '}';
    }
}
