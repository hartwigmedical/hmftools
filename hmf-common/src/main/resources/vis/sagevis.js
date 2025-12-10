let VISIBLE_READ_INFO = null;
let VISIBLE_FRAG_READS = null;

function main()
{
    $(".read-svg").on("click", function () {
        let readInfo = $($(this).siblings(".read-info"))[0];
        let fragmentReads = $(this).siblings(".read-of-fragment-sgv");
        if (VISIBLE_READ_INFO === null)
        {
            VISIBLE_READ_INFO = readInfo;
            VISIBLE_FRAG_READS = fragmentReads;
            $(readInfo).css("display", "block");
            $(fragmentReads).css("display", "block");
        }
        else if (readInfo.isSameNode(VISIBLE_READ_INFO))
        {
            $(VISIBLE_READ_INFO).css("display", "none");
            $(VISIBLE_FRAG_READS).css("display", "none");
            VISIBLE_READ_INFO = null;
            VISIBLE_FRAG_READS = null;
        }
        else
        {
            $(VISIBLE_READ_INFO).css("display", "none");
            $(VISIBLE_FRAG_READS).css("display", "none");
            VISIBLE_READ_INFO = readInfo;
            VISIBLE_FRAG_READS = fragmentReads;
            $(readInfo).css("display", "block");
            $(fragmentReads).css("display", "block");
        }
    });
}

main();