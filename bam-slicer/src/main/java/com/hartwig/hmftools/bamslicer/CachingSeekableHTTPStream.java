package com.hartwig.hmftools.bamslicer;

import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Chunk;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.HttpUtils;

public class CachingSeekableHTTPStream extends SeekableStream {
    private static final Logger LOGGER = LogManager.getLogger(CachingSeekableHTTPStream.class);
    private byte[] currentBytes = null;
    private long currentBytesOffset = 0;
    private long position = 0;
    private long contentLength = -1;
    private final ChunkHttpBuffer chunkBuffer;

    CachingSeekableHTTPStream(final URL url, @NotNull final List<Chunk> chunks, final int maxBufferSize) throws IOException {
        // Try to get the file length
        // Note: This also sets setDefaultUseCaches(false), which is important
        final String contentLengthString = HttpUtils.getHeaderField(url, "Content-Length");
        if (contentLengthString != null) {
            try {
                contentLength = Long.parseLong(contentLengthString);
            } catch (NumberFormatException ignored) {
                System.err.println("WARNING: Invalid content length (" + contentLengthString + "  for: " + url);
                contentLength = -1;
            }
        }
        LOGGER.info("Caching bam chunks from {}", url);
        chunkBuffer = new ChunkHttpBuffer(url, maxBufferSize, chunks);
        updatePosition(0);
    }

    @Override
    public long position() {
        return position;
    }

    @Override
    public long length() {
        return contentLength;
    }

    @Override
    public long skip(long n) throws IOException {
        long bytesToSkip = Math.min(n, contentLength - position);
        updatePosition(position + bytesToSkip);
        return bytesToSkip;
    }

    private void updatePosition(final long position) throws IOException {
        if (currentBytes == null || position >= currentBytesOffset + currentBytes.length) {
            final Map.Entry<Long, byte[]> bytesEntry = chunkBuffer.getEntryAtPosition(position);
            currentBytesOffset = bytesEntry.getKey();
            currentBytes = bytesEntry.getValue();
        }
        this.position = position;
    }

    @Override
    public boolean eof() throws IOException {
        return contentLength > 0 && position >= contentLength;
    }

    @Override
    public void seek(final long position) throws IOException {
        updatePosition(position);
    }

    @Override
    public int read(byte[] buffer, int offset, int len) throws IOException {
        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException("Offset=" + offset + ",len=" + len + ",buflen=" + buffer.length);
        }
        if (len == 0 || position == contentLength) {
            return 0;
        }
        final int sourcePosition = (int) (position - currentBytesOffset);
        System.arraycopy(currentBytes, sourcePosition, buffer, offset, len);
        updatePosition(position + len);
        return len;
    }

    @Override
    public void close() throws IOException {
        chunkBuffer.closeHttpClient();
    }

    @Override
    public int read() throws IOException {
        byte[] tmp = new byte[1];
        read(tmp, 0, 1);
        return (int) tmp[0] & 0xFF;
    }

    @Override
    public String getSource() {
        return chunkBuffer.url().toString();
    }
}
