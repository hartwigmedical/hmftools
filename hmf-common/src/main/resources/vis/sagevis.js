let VISIBLE_READ_INFO = null;
let VISIBLE_FRAG_READS = null;

function setDisplay(elements, display)
{
    elements.forEach(function(el) { el.style.display = display; });
}

function main()
{
    document.querySelectorAll(".read-svg").forEach(function(el) {
        el.addEventListener("click", function() {
            let siblings = Array.from(this.parentElement.children).filter(function(child) { return child !== this; }, this);
            let readInfo = siblings.find(function(child) { return child.matches(".read-info"); });
            let fragmentReads = siblings.filter(function(child) { return child.matches(".read-of-fragment-sgv"); });
            if (VISIBLE_READ_INFO === null)
            {
                VISIBLE_READ_INFO = readInfo;
                VISIBLE_FRAG_READS = fragmentReads;
                readInfo.style.display = "block";
                setDisplay(fragmentReads, "block");
            }
            else if (readInfo.isSameNode(VISIBLE_READ_INFO))
            {
                VISIBLE_READ_INFO.style.display = "none";
                setDisplay(VISIBLE_FRAG_READS, "none");
                VISIBLE_READ_INFO = null;
                VISIBLE_FRAG_READS = null;
            }
            else
            {
                VISIBLE_READ_INFO.style.display = "none";
                setDisplay(VISIBLE_FRAG_READS, "none");
                VISIBLE_READ_INFO = readInfo;
                VISIBLE_FRAG_READS = fragmentReads;
                readInfo.style.display = "block";
                setDisplay(fragmentReads, "block");
            }
        });
    });
}

main();