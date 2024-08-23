import { Elysia, error, t } from "elysia"
import { createClient } from "redis"
import { cors } from '@elysiajs/cors'

const redis = createClient()

redis.connect()

const app = new Elysia()
  .use(cors())
  .get("/api/result", async({ set, headers }) => {
    set.headers["access-control-allow-origin"] = "*"

    console.log(headers)
    const lastUpdate = await redis.get("us:last_update_at")
    if (lastUpdate === null) return error(404, "Not Found")
    else if(parseFloat(lastUpdate) <= headers.last_update_at) return error(304, "Not Modified")

    const result = await redis.lRange("us", 0, -1)
    return {
      lastUpdateAt: parseFloat(lastUpdate),
      result: result,
    }
  }, {
    headers: t.Object({
      last_update_at: t.Number(),
    })
  })
  .listen(3000)

console.log(
  `🦊 Elysia is running at ${app.server?.hostname}:${app.server?.port}`
);

export type App = typeof app