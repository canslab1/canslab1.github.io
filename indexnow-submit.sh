#!/usr/bin/env bash
# ──────────────────────────────────────────────
#  IndexNow 批次提交腳本
#  用途：通知搜尋引擎（Bing、Yandex、Naver 等）網站內容已更新
#  使用方式：
#    提交所有頁面：  ./indexnow-submit.sh
#    提交單一頁面：  ./indexnow-submit.sh https://canslab1.github.io/stories.html
# ──────────────────────────────────────────────

HOST="canslab1.github.io"
KEY="d22a81b36ccb45e085fe6679a822df52"
KEY_LOCATION="https://${HOST}/${KEY}.txt"
ENDPOINT="https://api.indexnow.org/IndexNow"

# 預設提交的頁面清單
ALL_URLS=(
    "https://${HOST}/"
    "https://${HOST}/index.html"
    "https://${HOST}/stories.html"
    "https://${HOST}/llms.txt"
)

# 如果有傳入參數，只提交指定的 URL
if [ $# -gt 0 ]; then
    URLS=("$@")
else
    URLS=("${ALL_URLS[@]}")
fi

echo "═══════════════════════════════════════"
echo "  IndexNow 批次提交"
echo "  Host: ${HOST}"
echo "  提交頁數: ${#URLS[@]}"
echo "═══════════════════════════════════════"

# 組合 URL 清單為 JSON 陣列
URL_JSON="["
for i in "${!URLS[@]}"; do
    if [ $i -gt 0 ]; then
        URL_JSON+=","
    fi
    URL_JSON+="\"${URLS[$i]}\""
done
URL_JSON+="]"

# 發送 POST 請求
RESPONSE=$(curl -s -o /dev/null -w "%{http_code}" -X POST "${ENDPOINT}" \
    -H "Content-Type: application/json; charset=utf-8" \
    -d "{
        \"host\": \"${HOST}\",
        \"key\": \"${KEY}\",
        \"keyLocation\": \"${KEY_LOCATION}\",
        \"urlList\": ${URL_JSON}
    }")

echo ""
echo "提交的 URL："
for url in "${URLS[@]}"; do
    echo "  ✓ ${url}"
done

echo ""
if [ "$RESPONSE" = "200" ]; then
    echo "✅ 提交成功！(HTTP ${RESPONSE})"
    echo "   搜尋引擎已收到更新通知。"
elif [ "$RESPONSE" = "202" ]; then
    echo "✅ 已接受！(HTTP ${RESPONSE})"
    echo "   搜尋引擎已接受提交，稍後處理。"
else
    echo "⚠️  回應碼：HTTP ${RESPONSE}"
    echo "   200=成功 | 202=已接受 | 422=無效URL | 429=頻率限制"
fi

echo ""
echo "═══════════════════════════════════════"
