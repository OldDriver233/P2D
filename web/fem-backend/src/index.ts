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
    const lastUpdate = await redis.get("last_update_at")
    if (lastUpdate === null) return error(404, "Not Found")
    else if(parseFloat(lastUpdate) <= headers.last_update_at) return error(304, "Not Modified")

    const u = await redis.lRange("u", 0, -1)
    const voltage = await redis.lRange("voltage", 0, -1)
    const delta_u = await redis.lRange("delta_u", 0, -1)
    const c_star = await redis.lRange("c_star", 0, -1)
    return {
      lastUpdateAt: parseFloat(lastUpdate),
      result: {
        u,
        voltage,
        delta_u,
        c_star,
      }
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