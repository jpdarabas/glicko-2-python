import math
import pandas as pd

def single_match_rating(player, score, opponent, tau = 0.2, rec = True):
  """
  Glicko-2 algorithm implementation in single matches. Based on Professor Mark E. Glickman's Glicko-2 system.

    Args:
        player: player/team DataFrame
        score: 0 for lose, 0.5 for draw and 1 for win
        opponent: opponent DataFrame
        tau: float constant

  """
  if player["ratings"] == None:
    player["ratings"] = {
        'rating': 1500,
        'RD': 350,
        'sigma': 0.06
    }
  player["g_ratings"] = {
      'mi' : (player.ratings['rating'] - 1500)/173.7178,
      'fi' : player.ratings['RD']/173.7178,
      'sigma' : player.ratings['sigma']
  }
  if opponent.ratings == None:
    opponent["ratings"] = {
        'rating': 1500,
        'RD': 350,
        'sigma': 0.06
    }
  opponent["g_ratings"] = {
      'mi' : (opponent.ratings['rating'] - 1500)/173.7178,
      'fi' : opponent.ratings['RD']/173.7178,
      'sigma' : opponent.ratings['sigma']
  }
  fi = player.g_ratings['fi']
  mi = player.g_ratings['mi']
  sigma = player.g_ratings['sigma']
  def g(fi):
    return 1/math.sqrt(1+3*(fi**2)/math.pi**2)
  def E(mij, fij, mi = player.g_ratings['mi']):
    return 1/(1+math.exp(-g(fij)*(mi-mij)))
  def quantity(score, player, opponent):
    mij = opponent.g_ratings['mi']
    fij = opponent.g_ratings['fi']
    v = g(fij)**2 * E(mij, fij) * (1 - E(mij, fij))
    v **= -1
    delta = g(opponent.g_ratings['fi'])*(score - E(opponent.g_ratings['mi'], opponent.g_ratings['fi']))
    delta *= v
    return v, delta
  v, delta = quantity(score, player, opponent)
  a = math.log(player.g_ratings['sigma']**2)
  def f(x):
    return ((math.e**x *(delta**2 - fi**2 - v - math.e**x))/((2*(fi**2 + v + math.e**x))**2) - (x-a)/(tau**2))
  y = 0.000001
  A = a
  if delta**2 > fi**2 + v:
    B = math.log(delta**2 - fi**2 - v)
  else:
    k = 1
    while f(a - k*tau)< 0:
      k += 1
    B = a - k*tau
  fa, fb = f(A), f(B)
  while abs(B - A) > y:
    C = A + (A-B)*fa/(fb-fa)
    fc = f(C)
    if fc*fb <= 0:
      A, fa = B, fb
    else:
      fa = fa/2
    B, fb = C, fc
  sigmaL = math.e**(A/2)

  fio = math.sqrt(fi**2 + sigmaL**2)
  fiL = 1/math.sqrt(1/(fio**2) + 1/v)
  
  miL = mi + fiL**2 * g(opponent.g_ratings['fi'])*(score - E(opponent.g_ratings['mi'], opponent.g_ratings['fi']))

  
  inverted_score = 1 if score == 0 else 0 if score == 1 else 0.5
  if rec == True:
    single_match_rating(opponent, inverted_score, player, tau, False)

  player["ratings"]['rating'] = 173.7178 * miL + 1500
  player["ratings"]['RD'] = 173.7178 * fiL
  player["ratings"]['sigma'] = sigmaL



## Testing

teams = [
    { 'ratings':{
        'rating':1500,
        'RD': 200,
        'sigma': 0.06
    },
        'g_ratings':None
    },
    { 'ratings':{
        'rating':1400,
        'RD': 30,
        'sigma': 0.06
    },
        'g_ratings':None
    },
    { 'ratings':{
        'rating':1550,
        'RD': 100,
        'sigma': 0.06
        },
        'g_ratings':None
    },
    { 'ratings':{
        'rating':1700,
        'RD': 300,
        'sigma': 0.06
        },
        'g_ratings':None
    },
]

teams_df = pd.DataFrame(teams)
single_match_rating(teams_df.iloc[0], 1, teams_df.iloc[1], 0.5)
single_match_rating(teams_df.iloc[0], 0, teams_df.iloc[2], 0.5)
single_match_rating(teams_df.iloc[0], 0, teams_df.iloc[3], 0.5)

print(teams_df.iloc[0]['ratings'])
# The results are approximately the same as the original algorithm, but not equals