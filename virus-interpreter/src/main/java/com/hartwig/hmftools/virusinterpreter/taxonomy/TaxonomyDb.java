package com.hartwig.hmftools.virusinterpreter.taxonomy;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class TaxonomyDb {

    private static final Logger LOGGER = LogManager.getLogger(TaxonomyDb.class);

    @NotNull
    private final Map<Integer, String> taxidToNameMap;

    public TaxonomyDb(@NotNull final Map<Integer, String> taxidToNameMap) {
        this.taxidToNameMap = taxidToNameMap;
    }

    @NotNull
    public String lookupName(int taxid) {
        boolean isMappedTaxid = taxidExists(taxid);

        if (!isMappedTaxid) {
            LOGGER.warn("Could not match taxid '{}' to a name in taxonomy DB", taxid);
        }

        return isMappedTaxid ? taxidToNameMap.get(taxid) : Strings.EMPTY;
    }

    @VisibleForTesting
    boolean taxidExists(int taxid) {
        return taxidToNameMap.containsKey(taxid);
    }

    @VisibleForTesting
    int count() {
        return taxidToNameMap.keySet().size();
    }

}
