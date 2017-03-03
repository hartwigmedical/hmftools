package com.hartwig.hmftools.fastqstats;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class FastqReader {
    private BufferedReader reader;
    private int offset = 33;
    private FastqTracker tracker;

    public FastqReader(InputStream in, FastqTracker tracker) {
        this.tracker = tracker;
        reader = new BufferedReader(new InputStreamReader(in));
    }

    private boolean hasNext() throws IOException {
        reader.mark(1);
        int c = reader.read();
        if (c == -1) {
            return false;
        }
        reader.reset();
        return true;
    }

    /**
     * Reads Fastq data and increments counters for the keys passed in as parameter
     *
     * @param keys Keys to track for this file
     * @throws IOException
     */
    public void read(TrackerKey[] keys) throws IOException {
        while (hasNext()) {
            // id
            reader.readLine();
            String seq = reader.readLine();
            // + line
            reader.readLine();
            readQualLine(seq.length(), keys);
        }
    }

    public void close() throws IOException {
        reader.close();
    }

    private void readQualLine(int expectedLength, TrackerKey[] keys) throws IOException {
        int numRead = 0;
        int c;
        while ((c = reader.read()) != -1 && c != '\r' && c != '\n') {
            numRead++;
            tracker.addValue(keys, c - offset);
        }
        if (numRead != expectedLength)
            throw new IOException("Length mismatch between quality and sequence lines.");
        handleCRLF(c);
    }

    private void handleCRLF(int c) throws IOException {
        if (c == '\r') {
            reader.mark(1);
            c = reader.read();
            if (c != '\n' && c != -1) {
                reader.reset();
            }
        }
    }
}
