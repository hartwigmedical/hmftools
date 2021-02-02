package com.hartwig.hmftools.patientdb.dao;

import java.sql.Timestamp;
import java.util.List;

interface BufferedWriterConsumer<T> {

    void initialise();

    void accept(Timestamp timestamp, List<T> entries);
}
