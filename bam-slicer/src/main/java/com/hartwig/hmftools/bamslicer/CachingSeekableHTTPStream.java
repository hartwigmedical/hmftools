package com.hartwig.hmftools.bamslicer;

import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Chunk;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import htsjdk.samtools.util.HttpUtils;
import io.reactivex.Observable;
import io.reactivex.subjects.AsyncSubject;
import io.reactivex.subjects.Subject;
import okhttp3.Call;
import okhttp3.Callback;
import okhttp3.ConnectionPool;
import okhttp3.Dispatcher;
import okhttp3.Headers;
import okhttp3.OkHttpClient;
import okhttp3.Request;
import okhttp3.Response;
import okhttp3.ResponseBody;

public class CachingSeekableHTTPStream extends SeekableStream {
    private static final Logger LOGGER = LogManager.getLogger(CachingSeekableHTTPStream.class);
    private static final int MAX_REQUESTS = 50;
    private final ConcurrentSkipListMap<Long, byte[]> bytesPerOffset = new ConcurrentSkipListMap<>();
    private byte[] currentBytes = null;
    private long currentBytesOffset = 0;
    private long position = 0;
    private long contentLength = -1;
    private final URL url;
    private final static OkHttpClient httpClient;

    static {
        final Dispatcher requestDispatcher = new Dispatcher();
        requestDispatcher.setMaxRequests(MAX_REQUESTS);
        requestDispatcher.setMaxRequestsPerHost(MAX_REQUESTS);
        httpClient = new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 1, TimeUnit.MINUTES))
                .dispatcher(requestDispatcher)
                .build();
    }

    CachingSeekableHTTPStream(final URL url, @NotNull final List<Chunk> chunks) {
        this.url = url;
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
        final Observable<Pair<Long, byte[]>> chunksBytes = Observable.fromIterable(chunks).flatMap(this::getBytesForChunk);
        chunksBytes.blockingSubscribe(pair -> bytesPerOffset.put(pair.getLeft(), pair.getRight()), LOGGER::error);
        LOGGER.info("Done caching...");
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

    private void updatePosition(final long position) {
        if (currentBytes == null || position >= currentBytesOffset + currentBytes.length) {
            final Map.Entry<Long, byte[]> bytesEntry = bytesPerOffset.floorEntry(position);
            currentBytesOffset = bytesEntry.getKey();
            currentBytes = bytesEntry.getValue();
            this.position = position;
        } else {
            this.position = position;
        }
    }

    @Override
    public boolean eof() throws IOException {
        return contentLength > 0 && position >= contentLength;
    }

    @Override
    public void seek(final long position) {
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
        httpClient.dispatcher().executorService().shutdown();
    }

    @Override
    public int read() throws IOException {
        byte[] tmp = new byte[1];
        read(tmp, 0, 1);
        return (int) tmp[0] & 0xFF;
    }

    @Override
    public String getSource() {
        return url.toString();
    }

    private Observable<Pair<Long, byte[]>> getBytesForChunk(@NotNull final Chunk chunk) {
        final long chunkStart = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkStart());
        final long chunkEnd = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkEnd());
        return getBytesBetweenOffsets(chunkStart, chunkEnd);
    }

    private Observable<Pair<Long, byte[]>> getBytesBetweenOffsets(final long start, final long end) {
        return Observable.defer(() -> {
            if (start <= end) {
                return readUrlBytes(start, end - start).map(bytes -> Pair.of(start, bytes));
            } else {
                return Observable.empty();
            }
        });
    }

    private Observable<byte[]> readUrlBytes(final long offset, final long count) {
        return Observable.defer(() -> {
            final Headers httpHeaders = new Headers.Builder().add("Range", "bytes=" + offset + "-" + (offset + count - 1)).build();
            final Request request = new Request.Builder().url(url).headers(httpHeaders).build();
            final Subject<byte[]> subject = AsyncSubject.create();
            final Callback callback = new Callback() {
                @Override
                public void onFailure(@NotNull final Call call, @NotNull final IOException e) {
                    subject.onError(e);
                }

                @Override
                public void onResponse(@NotNull final Call call, @NotNull final Response response) throws IOException {
                    final ResponseBody body = response.body();
                    try {
                        if (response.isSuccessful() && body != null) {
                            subject.onNext(body.bytes());
                            subject.onComplete();
                        } else {
                            subject.onError(new IOException("Response was not successful or body was null."));
                        }
                    } catch (final Exception e) {
                        subject.onError(e);
                    } finally {
                        if (body != null) {
                            body.close();
                        }
                    }
                }
            };
            httpClient.newCall(request).enqueue(callback);
            return subject;
        }).retry(5);
    }
}
