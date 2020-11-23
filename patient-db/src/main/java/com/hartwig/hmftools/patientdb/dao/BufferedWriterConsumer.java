package com.hartwig.hmftools.patientdb.dao;

import java.sql.Timestamp;
import java.util.List;

interface BufferedWriterConsumer<T> {
    void intialise();

    void accept(Timestamp timestamp, List<T> entries);
}
